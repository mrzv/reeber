add_executable              (tmt-mst        tmt-mst.cpp ${DEBUG_SOURCES} ${SOURCES})
target_link_libraries       (tmt-mst        ${libraries})

add_executable              (kruskal-mst    kruskal-mst.cpp ${DEBUG_SOURCES} ${SOURCES})
target_link_libraries       (kruskal-mst    ${libraries})
