cmake_minimum_required(VERSION 3.10)
project(BeamFDM LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)

add_executable(BeamFDM
    src/main.cpp
    src/BeamSolver.cpp
)

target_include_directories(BeamFDM PRIVATE ${EIGEN3_INCLUDE_DIR})