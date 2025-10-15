cmake_minimum_required(VERSION 3.10)
project(BeamFDM LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(BeamFDM
    src/main.cpp
    src/BeamSolver.cpp
)