cmake_minimum_required(VERSION 3.29.3)

project(
  "PoissonRecon"
  VERSION 1.0
  DESCRIPTION "An attempt at implementing Poisson reconstruction in C++"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS True)

include(FetchContent)

FetchContent_Declare(
  Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0)

FetchContent_MakeAvailable(Eigen)

add_executable(main src/main.cpp src/Octree.cpp src/utils.cpp src/io.cpp)
target_link_libraries(main PRIVATE Eigen3::Eigen)
add_subdirectory(tests)
