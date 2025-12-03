# ============================
#  COMPILERS & FLAGS
# ============================

# CPU build (OpenMP)
CXX = g++-15
CXXFLAGS = -O2 -std=c++17 -Wall -fopenmp

# GPU build (CUDA)
NVCC = nvcc
NVCCFLAGS = -O2 -std=c++17 -arch=sm_61

# ============================
#  TARGETS
# ============================

CPU_TARGET = nbody_omp
GPU_TARGET = nbody_cuda

CPU_SRC = src/main.cpp
GPU_SRC = src/cudaNBody.cpp

.PHONY: all omp cuda clean run_omp run_cuda srun_cuda

# Default: build everything
all: omp cuda

# ============================
#  CPU (OpenMP) BUILD
# ============================

omp: $(CPU_TARGET)

$(CPU_TARGET): $(CPU_SRC)
	$(CXX) $(CXXFLAGS) $(CPU_SRC) -o $(CPU_TARGET)

# ============================
#  GPU (CUDA) BUILD
# ============================

cuda: $(GPU_TARGET)

$(GPU_TARGET): $(GPU_SRC)
	$(NVCC) $(NVCCFLAGS) $(GPU_SRC) -o $(GPU_TARGET)

# ============================
#  RUN TARGETS
# ============================

# --- Run CPU version locally ---
run_omp: $(CPU_TARGET)
	./$(CPU_TARGET) sem 200 5000000 10000 > solar_omp.tsv
	python3 plot.py solar_omp.tsv solar_omp.pdf 10000

# --- Run CUDA version locally (ONLY if machine has a GPU) ---
run_cuda: $(GPU_TARGET)
	./$(GPU_TARGET) sem 200 5000000 10000 > solar_cuda.tsv
	python3 plot.py solar_cuda.tsv solar_cuda.pdf 10000

# --- Run CUDA version through SLURM on Centaurus ---
srun_cuda: $(GPU_TARGET)
	srun --partition=GPU --gres=gpu:1 ./$(GPU_TARGET) sem 200 5000000 10000 > solar_cuda.tsv
	python3 plot.py solar_cuda.tsv solar_cuda.pdf 10000


clean:
	rm -f $(CPU_TARGET) $(GPU_TARGET) *.o *.tsv *.pdf
