cmake_minimum_required(VERSION 3.31.0)

project(
  "PoissonRecon"
  VERSION 1.0
  DESCRIPTION "An attempt at implementing Poisson reconstruction in C++"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS True)
set(CMAKE_BUILD_TYPE Release)
# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_COMPILER "/opt/homebrew/opt/llvm/bin/clang")
set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
  set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=address")
endif()

cmake_policy(SET CMP0167 NEW)

# Add OpenMP include and library paths
set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include"
)
set(CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} -L/opt/homebrew/opt/libomp/lib -lomp")
set(OpenMP_CXX_FLAGS
    "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")

set(EIGEN_DIR "${CMAKE_SOURCE_DIR}/third_party/eigen")
include_directories(${EIGEN_DIR})

find_package(VTK REQUIRED)
find_package(GSL REQUIRED)

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(benchmarks)
