file(GLOB UTILS "utils/*.cpp")
list(REMOVE_ITEM UTILS ${CMAKE_CURRENT_SOURCE_DIR}/utils/plot.cpp)
set(PLOT ${CMAKE_CURRENT_SOURCE_DIR}/utils/plot.cpp)

file(GLOB NORMAL ${CMAKE_CURRENT_SOURCE_DIR}/Normal.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/Emst.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/RiemannianGraph.cpp)

add_executable(test test.cpp ${UTILS})
add_executable(
  main
  main.cpp
  PoissonRecon.cpp
  ${NORMAL}
  Octree.cpp
  pOctree.cpp
  HRefine.cpp
  ${UTILS})

target_link_libraries(main PUBLIC ${GSL_LIBRARIES})
add_executable(plotIndicator plotIndicator.cpp ${UTILS} ${PLOT})
target_link_libraries(plotIndicator PRIVATE ${VTK_LIBRARIES})

add_executable(plotScalar plotScalarField.cpp ${UTILS} ${PLOT})
target_link_libraries(plotScalar PRIVATE ${VTK_LIBRARIES})

add_executable(plotDebug Debug.cpp ${UTILS} ${PLOT})
target_link_libraries(plotDebug PRIVATE ${VTK_LIBRARIES})

# add_executable( plot_poisson plotPoisson.cpp PoissonRecon.cpp Normal.cpp
# RiemannianGraph.cpp Emst.cpp Octree.cpp pOctree.cpp ${UTILS})
# target_link_libraries(plot_poisson PRIVATE ${VTK_LIBRARIES} ${GSL_LIBRARIES})
#
# add_executable(plot plot.cpp Octree.cpp Normal.cpp RiemannianGraph.cpp
# Emst.cpp ${UTILS})
#
# target_link_libraries(plot PRIVATE ${VTK_LIBRARIES})
#
# add_executable(test test.cpp) target_link_libraries(test PRIVATE
# ${VTK_LIBRARIES})
