cmake_minimum_required(VERSION 3.16)
project(nodal_dg_1d LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Force CMake to export compile_commands.json for Neovim/LSP
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Include the headers directory
include_directories(inc)

# Define executables and their source files
add_executable(module_tests src/module_tests.cpp)
add_executable(linear_advection src/linear_advection.cpp)
