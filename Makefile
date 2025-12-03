# CPU build (OpenMP)
CXX = g++
CXXFLAGS = -O2 -std=c++17 -Wall -fopenmp

# GPU build (CUDA)
NVCC = nvcc
NVCCFLAGS = -O2 -std=c++17 -arch=sm_61


CPU_TARGET = nbody_omp
GPU_TARGET = nbody_cuda
SEQ_TARGET = nbody_seq

CPU_SRC = src/ompNBody.cpp
GPU_SRC = src/cudaNBody.cu
SEQ_SRC = src/seqNBody.cpp

.PHONY: all omp cuda clean run_omp run_cuda srun_cuda


all: omp cuda seq

seq: $(SEQ_TARGET)

$(SEQ_TARGET): $(SEQ_SRC)
	$(CXX) $(CXXFLAGS) $(SEQ_SRC) -o $(SEQ_TARGET)


omp: $(CPU_TARGET)

$(CPU_TARGET): $(CPU_SRC)
	$(CXX) $(CXXFLAGS) $(CPU_SRC) -o $(CPU_TARGET)


cuda: $(GPU_TARGET)

$(GPU_TARGET): $(GPU_SRC)
	$(NVCC) $(NVCCFLAGS) $(GPU_SRC) -o $(GPU_TARGET)



run_omp: $(CPU_TARGET)
	./$(CPU_TARGET) sem 200 900000000 10000 > solar_omp.tsv
	python3 plot.py solar_omp.tsv solar_omp.pdf 10000


run_cuda: $(GPU_TARGET)
	./$(GPU_TARGET) sem 200 5000000 10000 > solar_cuda.tsv
	python3 plot.py solar_cuda.tsv solar_cuda.pdf 10000

srun_cuda: $(GPU_TARGET)
	srun --partition=GPU --gres=gpu:1 ./$(GPU_TARGET) sem 200 5000000 10000 > solar_cuda.tsv
	python3 plot.py solar_cuda.tsv solar_cuda.pdf 10000


clean:
	rm -f $(CPU_TARGET) $(GPU_TARGET) *.o *.tsv *.pdf
