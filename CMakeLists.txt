cmake_minimum_required(VERSION 3.31.0)

project(
  "PoissonRecon"
  VERSION 1.0
  DESCRIPTION "An attempt at implementing Poisson reconstruction in C++"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS True)
set(CMAKE_BUILD_TYPE "Release")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
  set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=address")
endif()

cmake_policy(SET CMP0167 NEW)

# find_package(Boost REQUIRED)
find_package(VTK REQUIRED)
find_package(GSL REQUIRED)

set(EIGEN_DIR "${CMAKE_SOURCE_DIR}/third_party/eigen")
include_directories(${EIGEN_DIR})

add_subdirectory(src)
add_subdirectory(tests)
