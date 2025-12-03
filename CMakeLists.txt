cmake_minimum_required(VERSION 4.0)
project(N_bodySimulation)

set(CMAKE_CXX_STANDARD 20)

add_executable(N_bodySimulation
        src/ompNBody.cpp
        tooling/omp_loop.hpp
        tooling/omp_tasking.hpp
        tooling/seq_loop.hpp
        tooling/seq_tasking.hpp
        src/cudaNBody.cpp
        src/cudaNBody.cpp
)
