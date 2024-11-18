cmake_minimum_required(VERSION 3.29.3)

project(
  "PoissonRecon"
  VERSION 1.0
  DESCRIPTION "An attempt at implementing Poisson reconstruction in C++"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS True)

add_subdirectory(tests)
