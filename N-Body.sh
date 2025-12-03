#!/bin/bash
#SBATCH --job-name=nbody_gpu
#SBATCH --partition=GPU
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --output=job_output.txt

module load cuda/12.4


make all


date
./nbody_cuda 1000 200 500000 1000 256 > solar_cuda.out
date

date
./nbody_omp 1000 200 500000 1000 > solar_omp.out
date

date
./nbody_seq 1000 200 500000 1000 > solar_seq.out
date