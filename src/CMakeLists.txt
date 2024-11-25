file(GLOB UTILS "utils/*.cpp")

add_executable(
  main
  main.cpp
  PoissonRecon.cpp
  Normal.cpp
  RiemannianGraph.cpp
  Emst.cpp
  Octree.cpp
  pOctree.cpp
  basis.cpp
  integration.cpp
  ${UTILS})

add_executable(plotIndicator plotIndicator.cpp ${UTILS})
target_link_libraries(plotIndicator PRIVATE ${VTK_LIBRARIES})

# add_executable(main main.cpp basis.cpp integration.cpp ${UTILS})
target_link_libraries(main ${GSL_LIBRARIES})

add_executable(
  plot_poisson
  plotPoisson.cpp
  PoissonRecon.cpp
  Normal.cpp
  RiemannianGraph.cpp
  Emst.cpp
  Octree.cpp
  pOctree.cpp
  basis.cpp
  integration.cpp
  ${UTILS})
target_link_libraries(plot_poisson PRIVATE ${VTK_LIBRARIES} ${GSL_LIBRARIES})
# target_link_libraries(plot_poisson ${GSL_LIBRARIES})

#
# add_executable(plot plot.cpp Octree.cpp Normal.cpp RiemannianGraph.cpp
# Emst.cpp ${UTILS})
#
# target_link_libraries(plot PRIVATE ${VTK_LIBRARIES})
#
# add_executable(test test.cpp) target_link_libraries(test PRIVATE
# ${VTK_LIBRARIES})
