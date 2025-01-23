file(GLOB UTILS ${CMAKE_SOURCE_DIR}/src/utils/*.cpp)
list(REMOVE_ITEM UTILS ${CMAKE_SOURCE_DIR}/src/utils/plot.cpp)

include(FetchContent)

FetchContent_Declare(
  benchmark
  GIT_REPOSITORY https://github.com/google/benchmark.git
  GIT_TAG main)

set(BENCHMARK_USE_BUNDLED_GTEST OFF)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
set(BENCHMARK_ENABLE_TESTING OFF)
FetchContent_MakeAvailable(benchmark)

add_executable(
  benchmark_octree
  benchmark_octree.cpp ../src/Octree.cpp ../src/RiemannianGraph.cpp
  ../src/Emst.cpp ../src/Normal.cpp ${UTILS})
target_link_libraries(benchmark_octree PRIVATE benchmark::benchmark pthread)

add_executable(
  benchmark_poisson
  benchmark_poisson.cpp
  ../src/Octree.cpp
  ../src/RiemannianGraph.cpp
  ../src/Emst.cpp
  ../src/Normal.cpp
  ../src/PoissonRecon.cpp
  ../src/HRefine.cpp
  ${UTILS})
target_link_libraries(benchmark_poisson PRIVATE benchmark::benchmark pthread)
