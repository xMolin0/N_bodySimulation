## N-Body Simulation (CPU, OpenMP, CUDA)
This project implements an N-body gravitational simulation in C++ with three different execution models to explore parallel performance and scalability:
Sequential CPU implementation
Multi-threaded CPU implementation using OpenMP
GPU-accelerated implementation using CUDA
The goal of the project is to study how different hardware architectures and parallelization strategies impact the performance of an O(N²) physics simulation.
The simulator computes gravitational forces between bodies (planets, moons, stars, or randomly generated particles) and integrates their motion forward in time.

## Features
- Particle state includes: **mass, position, velocity, force**.
- Built-in initialization modes:
  - `sem` → Sun–Earth–Moon toy system
  - Random bodies (user chooses `N`)
- Gravitational force calculation with softening factor to prevent singularities.
- Euler integration of motion (`v += a*dt; x += v*dt`).
- Output in `.tsv` format for plotting with Python.
- Centering & scaling output (optional) so solar system fits nicely on charts.


## 🛠️ Build Instructions

### Prerequisites
- C++17 or newer (e.g., `g++`, `clang++`)
- `make` (for using the Makefile)
- Python 3 with `matplotlib` (for plotting)

### Build
Clone and build:
```bash
git clone https://github.com/xunit0/N-bodySimulation.git
cd N-bodySimulation
make run
```

### Example Workflow
```bash
make omp
make run_omp
open solar_omp.pdf
```
This builds the OpenMP version, runs the simulation, and generates the final visualization.

*Parameters for the executable can be changed in the makefile.
