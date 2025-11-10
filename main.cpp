// nbody.cpp  — simple student-style N-body
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <random>
#include <iomanip>


#include "tooling/omp_loop.hpp"

using namespace std;

struct Vec3 {
    double x=0, y=0, z=0;
};

struct Particle {
    double m=1.0;
    Vec3 pos, vel, force;
};

struct Config {
    double G = 6.674e-11;   // gravitational constant
    double dt = 1.0;        // time step
    int steps = 1000;       // number of iterations
    int dump_every = 10;    // how often to print a line
    double softening = 1e-3; // small value to avoid div-by-0 (unit-dependent)
    unsigned seed = 42;     // RNG for random init
};

struct NBodyState {
    vector<Particle> p;

    size_t size() const { return p.size(); }
    void resize(size_t n) { p.assign(n, Particle{}); }

    // --- init: random ---
    void init_random(size_t N, double mass_min, double mass_max,
                     double pos_span, double vel_span, unsigned seed=42) {
        resize(N);
        mt19937 rng(seed);
        uniform_real_distribution<double> Um(mass_min, mass_max);
        uniform_real_distribution<double> U(-1.0, 1.0);
        for (auto &b : p) {
            b.m = Um(rng);
            b.pos = { U(rng)*pos_span, U(rng)*pos_span, U(rng)*pos_span };
            b.vel = { U(rng)*vel_span, U(rng)*vel_span, U(rng)*vel_span };
            b.force = {0,0,0};
        }
    }


    void init_sem() {
        enum Planets {SUN, MERCURY, VENUS, EARTH, MARS, JUPITER, SATURN, URANUS, NEPTUNE, MOON};
        resize(10);

        // Masses in kg
        p[SUN].m = 1.9891 * std::pow(10, 30);
        p[MERCURY].m = 3.285 * std::pow(10, 23);
        p[VENUS].m = 4.867 * std::pow(10, 24);
        p[EARTH].m = 5.972 * std::pow(10, 24);
        p[MARS].m = 6.39 * std::pow(10, 23);
        p[JUPITER].m = 1.898 * std::pow(10, 27);
        p[SATURN].m = 5.683 * std::pow(10, 26);
        p[URANUS].m = 8.681 * std::pow(10, 25);
        p[NEPTUNE].m = 1.024 * std::pow(10, 26);
        p[MOON].m = 7.342 * std::pow(10, 22);

        // Positions (in meters) and velocities (in m/s)
        double AU = 1.496 * std::pow(10, 11); // Astronomical Unit

        p[SUN].pos     = {0, 0, 0};
        p[MERCURY].pos = {0.39 * AU, 0, 0};
        p[VENUS].pos   = {0.72 * AU, 0, 0};
        p[EARTH].pos   = {1.0 * AU,  0, 0};
        p[MARS].pos    = {1.52 * AU, 0, 0};
        p[JUPITER].pos = {5.20 * AU, 0, 0};
        p[SATURN].pos  = {9.58 * AU, 0, 0};
        p[URANUS].pos  = {19.22 * AU, 0, 0};
        p[NEPTUNE].pos = {30.05 * AU, 0, 0};
        p[MOON].pos    = {1.0 * AU + 3.844e8, 0, 0}; // offset from Earth

        // Velocities (m/s) — roughly tangential orbital speeds
        p[SUN].vel     = {0, 0, 0};
        p[MERCURY].vel = {0, 47870, 0};
        p[VENUS].vel   = {0, 35020, 0};
        p[EARTH].vel   = {0, 29780, 0};
        p[MARS].vel    = {0, 24130, 0};
        p[JUPITER].vel = {0, 13070, 0};
        p[SATURN].vel  = {0, 9680,  0};
        p[URANUS].vel  = {0, 6800,  0};
        p[NEPTUNE].vel = {0, 5430,  0};
        p[MOON].vel    = {0, 29780 + 1022, 0}; // Moon relative to Earth

        // Initialize forces to zero
        for (auto &b : p)
            b.force = {0, 0, 0};
    }


    bool load_from_file(const std::string& path) {
        std::ifstream fin(path);
        if (!fin) return false;

        // read the first non-empty, non-comment line
        std::string line;
        while (std::getline(fin, line)) {
            // trim leading spaces
            std::size_t i = 0;
            while (i < line.size() && std::isspace((unsigned char)line[i])) ++i;
            if (i >= line.size() || line[i] == '#') continue; // skip blank or comment
            line.erase(0, i);
            break;
        }
        if (line.empty()) return false;

        std::istringstream iss(line);

        int N;
        if (!(iss >> N) || N <= 0) return false;
        resize((std::size_t)N);

        for (int i = 0; i < N; ++i) {
            Particle b;
            if (!(iss >> b.m
                      >> b.pos.x >> b.pos.y >> b.pos.z
                      >> b.vel.x >> b.vel.y >> b.vel.z
                      >> b.force.x >> b.force.y >> b.force.z)) {
                return false; // not enough numbers on the line
                      }
            p[i] = b;
        }
        return true;
    }

};

// --- math helpers ---
static inline Vec3 add(const Vec3&a,const Vec3&b){ return {a.x+b.x,a.y+b.y,a.z+b.z}; }
static inline Vec3 sub(const Vec3&a,const Vec3&b){ return {a.x-b.x,a.y-b.y,a.z-b.z}; }
static inline Vec3 mul(const Vec3&a,double s){ return {a.x*s,a.y*s,a.z*s}; }

static inline double norm2(const Vec3& v){ return v.x*v.x + v.y*v.y + v.z*v.z; }

void compute_forces(NBodyState& S, const Config& cfg, OmpLoop& omp) {
    // reset
    for (auto &b : S.p) b.force = {0,0,0};

    const double eps2 = cfg.softening * cfg.softening;
    int N = (int)S.size();

    // TLS
    struct ThreadLocal {
        std::vector<Vec3> local_forces;
    };

    omp.parfor<ThreadLocal>(
        0, N, 1,

        // BEFORE
        [N](ThreadLocal& tls) {
            tls.local_forces.resize(N, {0, 0, 0});
        },

        // LOOP BODY
        [&](size_t i, ThreadLocal& tls) {
            for (int j = i+1; j < N; j++) {
                Vec3 rij = sub(S.p[j].pos, S.p[i].pos);
                double r2 = norm2(rij) + eps2;
                double r = sqrt(r2);
                if (r == 0) continue;

                double Fmag = cfg.G * S.p[i].m * S.p[j].m / r2;
                Vec3 rhat = mul(rij, 1.0/r);
                Vec3 Fij = mul(rhat, Fmag);


                tls.local_forces[i] = add(tls.local_forces[i], Fij);
                tls.local_forces[j] = sub(tls.local_forces[j], Fij);
            }
        },

        // AFTER
        [&S](ThreadLocal& tls) {

            for (size_t i = 0; i < S.size(); i++) {
                S.p[i].force = add(S.p[i].force, tls.local_forces[i]);
            }
        }
    );
}

// --- integrate one step (Euler: v += a*dt; x += v*dt) ---
void step_euler(NBodyState& S, const Config& cfg) {
    for (auto &b : S.p) {
        Vec3 a = { b.force.x / b.m, b.force.y / b.m, b.force.z / b.m };
        b.vel = add(b.vel, mul(a, cfg.dt));
        b.pos = add(b.pos, mul(b.vel, cfg.dt)); // uses v_new per spec
    }
}

void write_state_tsv(const NBodyState& S) {

    std::cout << S.size() << '\t';

    for (size_t i=0; i<S.size(); ++i) {
        std::cout<<S.p[i].m<<'\t';
        std::cout<<S.p[i].pos.x<<'\t'<<S.p[i].pos.y<<'\t'<<S.p[i].pos.z<<'\t';
        std::cout<<S.p[i].vel.x<<'\t'<<S.p[i].vel.y<<'\t'<<S.p[i].vel.z<<'\t';
        std::cout<<S.p[i].force.x<<'\t'<<S.p[i].force.y<<'\t'<<S.p[i].force.z<<'\t';
    }
    cout << '\n';
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);

    OmpLoop omp;
    omp.setNbThread(11);

    string mode = argv[1];
    Config cfg;
    cfg.dt = stod(argv[2]);
    cfg.steps = stoi(argv[3]);
    cfg.dump_every = stoi(argv[4]);


    NBodyState S;

    auto is_number = [](const std::string& s) -> bool {
        if (s.empty()) return false;
        for (unsigned char c : s) {
            if (!std::isdigit(c)) return false;
        }
        return true;
    };


    if (mode == "sem") {
        S.init_sem();
    } else if (is_number(mode)) {
        size_t N = (size_t)stoull(mode);
        // simple random ranges: pick something sane (units arbitrary)
        S.init_random(N,
              1e22,   // mass_min: ~small dwarf planet
              1e27,   // mass_max: ~gas giant like Jupiter
              1.0e11, // pos_span: ~1 AU (distance Earth–Sun)
              1.0e4,  // vel_span: ~10 km/s (planet orbital speeds)
              cfg.seed);

    } else {
        if (!S.load_from_file(mode)) {
            cerr << "failed to load file: " << mode << "\n";
            return 2;
        }
    }


    compute_forces(S, cfg, omp);
    write_state_tsv(S);

    for (int step=1; step<=cfg.steps; ++step) {
        // 1) forces at current positions
        compute_forces(S, cfg, omp);
        // 2) integrate
        step_euler(S, cfg);
        // 3) output occasionally
        if (step % cfg.dump_every == 0) {
            // recompute forces for logging (forces correspond to printed state)
            compute_forces(S, cfg, omp);
            write_state_tsv(S);
        }
    }
    return 0;
}
