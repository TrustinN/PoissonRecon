file(GLOB UTILS "utils/*.cpp")

add_executable(
  main
  main.cpp
  Octree.cpp
  # Normal.cpp
  # RiemannianGraph.cpp
  # Emst.cpp
  # pOctree.cpp
  Node.cpp
  ${UTILS})

add_executable(
  plot
  plot.cpp
  Octree.cpp
  Normal.cpp
  RiemannianGraph.cpp
  Emst.cpp
  Node.cpp
  ${UTILS})

target_link_libraries(plot PRIVATE ${VTK_LIBRARIES})

add_executable(test test.cpp)
target_link_libraries(test PRIVATE ${VTK_LIBRARIES})
