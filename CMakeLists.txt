cmake_minimum_required(VERSION 3.29.3)

project(
  "PoissonRecon"
  VERSION 1.0
  DESCRIPTION "An attempt at implementing Poisson reconstruction in C++"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS True)

set(EIGEN_DIR "${CMAKE_SOURCE_DIR}/third_party/eigen")
include_directories(${EIGEN_DIR})

add_executable(
  main
  src/main.cpp
  src/Octree.cpp
  src/utils.cpp
  src/io.cpp
  src/sampling.cpp
  src/Normal.cpp
  src/RiemannianGraph.cpp
  src/Emst.cpp)
add_subdirectory(tests)
