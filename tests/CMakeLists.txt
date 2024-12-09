include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
FetchContent_MakeAvailable(googletest)

enable_testing()
include(GoogleTest)

file(GLOB UTILS "../src/utils/*.cpp")
list(REMOVE_ITEM UTILS ${CMAKE_CURRENT_SOURCE_DIR}/../src/utils/plot.cpp)
message(${UTILS})

add_executable(octree octree.cpp ../src/Octree.cpp ${UTILS})
target_link_libraries(octree GTest::gtest_main)
gtest_discover_tests(octree)

add_executable(riemannian riemannian.cpp ../src/Octree.cpp
                          ../src/RiemannianGraph.cpp ${UTILS})
target_link_libraries(riemannian GTest::gtest_main)
gtest_discover_tests(riemannian)

add_executable(normal normal.cpp ../src/Normal.cpp ../src/RiemannianGraph.cpp
                      ../src/Emst.cpp ../src/Octree.cpp ${UTILS})
target_link_libraries(normal GTest::gtest_main)
gtest_discover_tests(normal)

add_executable(polynomial polynomial.cpp ${UTILS})
target_link_libraries(polynomial GTest::gtest_main)
gtest_discover_tests(polynomial)

add_executable(ppolynomial ppolynomial.cpp ${UTILS})
target_link_libraries(ppolynomial GTest::gtest_main)
gtest_discover_tests(ppolynomial)

add_executable(bspline bspline.cpp ${UTILS})
target_link_libraries(bspline GTest::gtest_main)
gtest_discover_tests(bspline)

add_executable(polynomialxd polynomialxd.cpp ${UTILS})
target_link_libraries(polynomialxd GTest::gtest_main)
gtest_discover_tests(polynomialxd)
