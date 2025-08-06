# Specify the minimum version of CMake required.
cmake_minimum_required(VERSION 3.16)

# Define the project name.
project(Mie)

# Set the C++ standard to C++20. This is needed for the <numbers> header.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the executable target and specify its source files.
# This tells CMake to compile example.cpp and bhmie.cpp together
# to create an executable named "mie_example".
add_executable(mie_example example.cpp bhmie.cpp)


find_package(pybind11 REQUIRED)

pybind11_add_module(bhmie_py pybind_bhmie.cpp bhmie.cpp)
