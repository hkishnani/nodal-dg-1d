cmake_minimum_required(VERSION 3.16)
project(nodal_dg_1d LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Force CMake to export compile_commands.json for Neovim/LSP
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Include the headers directory
include_directories(inc)

# Define executables and their source files
add_executable(linear_advection src/linear_advection.cpp)

add_executable(module_test_mass_matrix src/module_test_mass_matrix.cpp)

add_executable(module_test_quadrature_weights_and_abscissae src/module_test_quadrature_weights_and_abscissae.cpp)

add_executable(module_test_basis_function src/module_test_basis_function.cpp)
