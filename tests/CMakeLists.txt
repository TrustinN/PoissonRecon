include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
FetchContent_MakeAvailable(googletest)

enable_testing()
include(GoogleTest)

file(GLOB UTILS "../src/utils/*.cpp")

add_executable(octree tests_octree.cpp ../src/Octree.cpp ${UTILS})
target_link_libraries(octree GTest::gtest_main ${VTK_LIBRARIES})
gtest_discover_tests(octree)

add_executable(
  normal_utils
  test_normal_utils.cpp ../src/Normal.cpp ../src/RiemannianGraph.cpp
  ../src/Emst.cpp ../src/Octree.cpp ${UTILS})
target_link_libraries(normal_utils GTest::gtest_main ${VTK_LIBRARIES})
gtest_discover_tests(normal_utils)
