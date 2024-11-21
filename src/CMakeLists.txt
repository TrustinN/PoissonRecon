add_executable(
  main
  main.cpp
  Octree.cpp
  utils/utils.cpp
  utils/io.cpp
  utils/sampling.cpp
  Normal.cpp
  RiemannianGraph.cpp
  Emst.cpp)

add_executable(
  plot
  plot.cpp
  Octree.cpp
  utils/utils.cpp
  utils/io.cpp
  utils/sampling.cpp
  Normal.cpp
  RiemannianGraph.cpp
  Emst.cpp)

target_link_libraries(plot PRIVATE ${VTK_LIBRARIES})

add_executable(test test.cpp)
target_link_libraries(test PRIVATE ${VTK_LIBRARIES})
