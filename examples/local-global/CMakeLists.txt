foreach                     (real   float double)

add_executable              (mt-np-local-global-${real}   mt-np-local-global.cpp ${DEBUG_SOURCES} ${SOURCES})
target_link_libraries       (mt-np-local-global-${real}   ${libraries})
set_target_properties       (mt-np-local-global-${real}   PROPERTIES COMPILE_FLAGS -DREEBER_REAL=${real})

endforeach                  (real)
