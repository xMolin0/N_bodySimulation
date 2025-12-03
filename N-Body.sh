#!/bin/bash
#SBATCH --job-name=nbody_gpu
#SBATCH --partition=GPU
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --output=job_output.txt

module load cuda/12.4

# Build the GPU version (optional)
make cuda

# === ACTUAL COMMAND TO EXECUTE ===
./nbody_cuda 10000 200 5000000 10000 > solar_cuda.out
