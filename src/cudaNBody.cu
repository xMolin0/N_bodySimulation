//
// Created by Xavier Molina on 12/2/25.
//

#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <cuda_runtime.h>

__device__ double G_dev = 6.674e-11;

struct simulation {
    size_t nbpart;

    double *d_mass, *d_x, *d_y, *d_z;
    double *d_vx, *d_vy, *d_vz;
    double *d_fx, *d_fy, *d_fz;

    // Constructor to allocate GPU memory
    simulation(size_t n) : nbpart(n) {
        size_t bytes = n * sizeof(double);
        cudaMalloc(&d_mass, bytes);
        cudaMalloc(&d_x, bytes);
        cudaMalloc(&d_y, bytes);
        cudaMalloc(&d_z, bytes);
        cudaMalloc(&d_vx, bytes);
        cudaMalloc(&d_vy, bytes);
        cudaMalloc(&d_vz, bytes);
        cudaMalloc(&d_fx, bytes);
        cudaMalloc(&d_fy, bytes);
        cudaMalloc(&d_fz, bytes);
    }

    // Destructor to free GPU memory
    ~simulation() {
        cudaFree(d_mass);
        cudaFree(d_x);
        cudaFree(d_y);
        cudaFree(d_z);
        cudaFree(d_vx);
        cudaFree(d_vy);
        cudaFree(d_vz);
        cudaFree(d_fx);
        cudaFree(d_fy);
        cudaFree(d_fz);
    }
};




__global__
void reset_forces(double* fx, double* fy, double* fz, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        fx[i] = 0;
        fy[i] = 0;
        fz[i] = 0;
    }
}

__global__
void compute_forces(double* mass, double* x, double* y, double* z,
                    double* fx, double* fy, double* fz,
                    int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    double xi = x[i];
    double yi = y[i];
    double zi = z[i];
    double mi = mass[i];

    double Fx = 0, Fy = 0, Fz = 0;

    for (int j = 0; j < n; j++) {
        if (i == j) continue;

        double dx = x[j] - xi;
        double dy = y[j] - yi;
        double dz = z[j] - zi;

        double distSq = dx*dx + dy*dy + dz*dz + 0.1;
        double dist = sqrt(distSq);

        double F = G_dev * mi * mass[j] / distSq;

        Fx += F * dx / dist;
        Fy += F * dy / dist;
        Fz += F * dz / dist;
    }

    fx[i] = Fx;
    fy[i] = Fy;
    fz[i] = Fz;
}

__global__
void integrate(double* x, double* y, double* z,
               double* vx, double* vy, double* vz,
               double* fx, double* fy, double* fz,
               double* mass,
               double dt,
               int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    vx[i] += fx[i] / mass[i] * dt;
    vy[i] += fy[i] / mass[i] * dt;
    vz[i] += fz[i] / mass[i] * dt;

    x[i] += vx[i] * dt;
    y[i] += vy[i] * dt;
    z[i] += vz[i] * dt;
}

void init_solar(simulation& s) {
    const int N = 10;
    s = simulation(N); // allocate GPU memory

    // CPU arrays for masses, positions, velocities
    double h_mass[N] = {
        1.9891e30, 3.285e23, 4.867e24, 5.972e24, 6.39e23,
        1.898e27, 5.683e26, 8.681e25, 1.024e26, 7.342e22
    };

    double AU = 1.496e11;

    double h_x[N] = {0, 0.39*AU, 0.72*AU, 1.0*AU, 1.52*AU, 5.20*AU, 9.58*AU, 19.22*AU, 30.05*AU, 1.0*AU + 3.844e8};
    double h_y[N] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double h_z[N] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    double h_vx[N] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double h_vy[N] = {0, 47870, 35020, 29780, 24130, 13070, 9680, 6800, 5430, 29780+1022};
    double h_vz[N] = {0,0,0,0,0,0,0,0,0,0};

    // Copy to GPU memory
    cudaMemcpy(s.d_mass, h_mass, N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(s.d_x, h_x, N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(s.d_y, h_y, N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(s.d_z, h_z, N*sizeof(double), cudaMemcpyHostToDevice);

    cudaMemcpy(s.d_vx, h_vx, N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(s.d_vy, h_vy, N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(s.d_vz, h_vz, N*sizeof(double), cudaMemcpyHostToDevice);
}

void dump_state(simulation& s, size_t N) {
    std::vector<double> mass(N), x(N), y(N), z(N);
    std::vector<double> vx(N), vy(N), vz(N);
    std::vector<double> fx(N), fy(N), fz(N);

    cudaMemcpy(mass.data(), s.d_mass, N*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(x.data(), s.d_x, N*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(y.data(), s.d_y, N*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(z.data(), s.d_z, N*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(vx.data(), s.d_vx, N*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(vy.data(), s.d_vy, N*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(vz.data(), s.d_vz, N*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(fx.data(), s.d_fx, N*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(fy.data(), s.d_fy, N*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(fz.data(), s.d_fz, N*sizeof(double), cudaMemcpyDeviceToHost);

    std::cout<<s.nbpart<<'\t';
    for (size_t i=0; i<N; ++i) {
        std::cout<<mass[i]<<'\t';
        std::cout<<x[i]<<'\t'<<y[i]<<'\t'<<z[i]<<'\t';
        std::cout<<vx[i]<<'\t'<<vy[i]<<'\t'<<vz[i]<<'\t';
        std::cout<<fx[i]<<'\t'<<fy[i]<<'\t'<<fz[i]<<'\t';
    }
    std::cout<<'\n';
}

void random_init(simulation& s) {
    size_t N = s.nbpart;
    // Allocate CPU arrays
    std::vector<double> mass(N), x(N), y(N), z(N);
    std::vector<double> vx(N), vy(N), vz(N);
    std::vector<double> fx(N), fy(N), fz(N);



    std::mt19937 rng(1234);                       // fixed seed for reproducibility
    std::normal_distribution<float> norm(0.0f, 1.0f);
    std::uniform_real_distribution<float> vel(-0.001f, 0.001f);

    for (int i = 0; i < N; i++) {
        mass[i] = 1.0f;

        // Gaussian distributed positions
        x[i] = norm(rng);
        y[i] = norm(rng);
        z[i] = norm(rng);

        // tiny random initial velocity
        vx[i] = vel(rng);
        vy[i] = vel(rng);
        vz[i] = vel(rng);
    }

    size_t bytes = N * sizeof(double);
    // Copy initial state
    cudaMemcpy(s.d_mass, mass.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(s.d_x, x.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(s.d_y, y.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(s.d_z, z.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(s.d_vx, vx.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(s.d_vy, vy.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(s.d_vz, vz.data(), bytes, cudaMemcpyHostToDevice);

}

int main(int argc, char* argv[]) {
    if (argc != 6) {
        std::cerr
          <<"usage: "<<argv[0]<<" <input> <dt> <nbstep> <printevery>"<<"\n"
          <<"input can be:"<<"\n"
          <<"a number (random initialization)"<<"\n"
          <<"planet (initialize with solar system)"<<"\n"
          <<"a filename (load from file in singleline tsv)"<<"\n";
        return -1;
    }

    double dt = std::atof(argv[2]); //in seconds
    size_t steps = std::atol(argv[3]);
    size_t printevery = std::atol(argv[4]);
    int block = std::atol(argv[5]);


    simulation s(1);
    //parse command line
    size_t N = std::atol(argv[1]); //return 0 if not a number
    if (N > 0) {
        s = simulation(N);
        random_init(s);
    } else {
        std::string inputparam = argv[1];
        if (inputparam == "planet") {
            init_solar(s);
        } else{
            //load_from_file(s, inputparam);
        }
    }


    int grid = (N + block - 1) / block;
    for (__uint32_t step = 0; step < steps; step++) {
        if (step % printevery == 0) dump_state(s, s.nbpart);

        reset_forces<<<grid, block>>>(s.d_fx, s.d_fy, s.d_fz, N);
        compute_forces<<<grid, block>>>(s.d_mass, s.d_x, s.d_y, s.d_z,
                                        s.d_fx, s.d_fy, s.d_fz, N);
        integrate<<<grid, block>>>(s.d_x, s.d_y, s.d_z,
                                   s.d_vx, s.d_vy, s.d_vz,
                                   s.d_fx, s.d_fy, s.d_fz,
                                   s.d_mass,
                                   dt, N);
    }



    return 0;
}
