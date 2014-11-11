foreach                     (real   float double)

add_executable              (mt-global-${real}   mt-global.cpp format.cc ${DEBUG_SOURCES})
target_link_libraries       (mt-global-${real}   ${libraries})
set_target_properties       (mt-global-${real}   PROPERTIES COMPILE_FLAGS -DREEBER_REAL=${real})

add_executable              (mt-local-global-${real}   mt-local-global.cpp format.cc ${DEBUG_SOURCES})
target_link_libraries       (mt-local-global-${real}   ${libraries})
set_target_properties       (mt-local-global-${real}   PROPERTIES COMPILE_FLAGS -DREEBER_REAL=${real})

add_executable              (persistence-${real}   persistence.cpp format.cc ${DEBUG_SOURCES})
target_link_libraries       (persistence-${real}   ${libraries})
set_target_properties       (persistence-${real}   PROPERTIES COMPILE_FLAGS -DREEBER_REAL=${real})

endforeach                  (real)
