foreach                     (real   float double)

add_executable              (mt-np-local-global-${real}   mt-np-local-global.cpp ${DEBUG_SOURCES} ${SOURCES})
target_link_libraries       (mt-np-local-global-${real}   ${libraries})
set_target_properties       (mt-np-local-global-${real}   PROPERTIES COMPILE_FLAGS -DREEBER_REAL=${real})

add_executable              (mt-np-lg-ghosts-${real}   mt-np-lg-ghosts.cpp ${DEBUG_SOURCES} ${SOURCES})
target_link_libraries       (mt-np-lg-ghosts-${real}   ${libraries})
set_target_properties       (mt-np-lg-ghosts-${real}   PROPERTIES COMPILE_FLAGS -DREEBER_REAL=${real})

add_executable              (persistence-lg-${real}   persistence-lg.cpp ${DEBUG_SOURCES} ${SOURCES})
target_link_libraries       (persistence-lg-${real}   ${libraries})
set_target_properties       (persistence-lg-${real}   PROPERTIES COMPILE_FLAGS -DREEBER_REAL=${real})

add_executable              (persistent-integral-lg-${real}   persistent-integral-lg.cpp ${DEBUG_SOURCES} ${SOURCES})
target_link_libraries       (persistent-integral-lg-${real}   ${libraries})
set_target_properties       (persistent-integral-lg-${real}   PROPERTIES COMPILE_FLAGS -DREEBER_REAL=${real})

endforeach                  (real)

if                          (boxlib)
string                      (TOLOWER ${CCSE_PRECISION} CCSE_PRECISION_LC)
add_executable              (mt-bl-local-global-${CCSE_PRECISION_LC}-${CCSE_SPACEDIM}d   mt-bl-local-global.cpp ${DEBUG_SOURCES} ${SOURCES})
target_link_libraries       (mt-bl-local-global-${CCSE_PRECISION_LC}-${CCSE_SPACEDIM}d   ${libraries} ${CCSE_LIBRARIES})
endif                       (boxlib)
