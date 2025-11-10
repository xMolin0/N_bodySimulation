# Makefile for N-Body Simulation

CXX = g++-15
CXXFLAGS = -O2 -std=c++17 -Wall -fopenmp

TARGET = nbody
SRC = main.cpp

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

run: $(TARGET)
	./$(TARGET) 3 100 100000 100 > solar.tsv
	python3 plot.py solar.tsv solar.pdf 10000

clean:
	rm -f $(TARGET) *.o *.tsv *.pdf
