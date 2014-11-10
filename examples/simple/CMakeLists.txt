add_executable              (simple-float   simple.cpp format.cc ${DEBUG_SOURCES})
target_link_libraries       (simple-float   ${libraries})
set_target_properties       (simple-float   PROPERTIES COMPILE_FLAGS -DREEBER_REAL=float)

add_executable              (simple-double  simple.cpp format.cc ${DEBUG_SOURCES})
target_link_libraries       (simple-double  ${libraries})
set_target_properties       (simple-double  PROPERTIES COMPILE_FLAGS -DREEBER_REAL=double)
