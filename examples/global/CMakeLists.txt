foreach                     (real   float double)

add_executable              (mt-np-global-${real}   mt-np-global.cpp ${DEBUG_SOURCES} ${SOURCES})
target_link_libraries       (mt-np-global-${real}   ${libraries})
set_target_properties       (mt-np-global-${real}   PROPERTIES COMPILE_FLAGS -DREEBER_REAL=${real})

add_executable              (tmt-np-global-${real}   tmt-np-global.cpp ${DEBUG_SOURCES} ${SOURCES})
target_link_libraries       (tmt-np-global-${real}   ${libraries})
set_target_properties       (tmt-np-global-${real}   PROPERTIES COMPILE_FLAGS -DREEBER_REAL=${real})

add_executable              (persistence-${real}   persistence.cpp ${DEBUG_SOURCES} ${SOURCES})
target_link_libraries       (persistence-${real}   ${libraries})
set_target_properties       (persistence-${real}   PROPERTIES COMPILE_FLAGS -DREEBER_REAL=${real})

endforeach                  (real)
