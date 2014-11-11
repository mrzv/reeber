foreach                     (real   float double)

add_executable              (simple-${real}   simple.cpp format.cc ${DEBUG_SOURCES})
target_link_libraries       (simple-${real}   ${libraries})
set_target_properties       (simple-${real}   PROPERTIES COMPILE_FLAGS -DREEBER_REAL=${real})

add_executable              (persistence-${real}   persistence.cpp format.cc ${DEBUG_SOURCES})
target_link_libraries       (persistence-${real}   ${libraries})
set_target_properties       (persistence-${real}   PROPERTIES COMPILE_FLAGS -DREEBER_REAL=${real})

endforeach                  (real)
