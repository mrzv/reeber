option                      (persistent-integral-trace-vertices  "Record vertices belonging to persistent integral components" OFF)

foreach                     (real   float double)

add_executable              (mt-local-global-${real}   mt-local-global.cpp ${DEBUG_SOURCES} ${SOURCES})
target_link_libraries       (mt-local-global-${real}   ${libraries})
set_target_properties       (mt-local-global-${real}   PROPERTIES COMPILE_FLAGS "-DREEBER_REAL=${real}")

add_executable              (mt-lg-ghosts-${real}     mt-lg-ghosts.cpp ${DEBUG_SOURCES} ${SOURCES})
target_link_libraries       (mt-lg-ghosts-${real}     ${libraries})
set_target_properties       (mt-lg-ghosts-${real}     PROPERTIES COMPILE_FLAGS "-DREEBER_REAL=${real}")

add_executable              (persistence-lg-${real}   persistence-lg.cpp ${DEBUG_SOURCES} ${SOURCES})
target_link_libraries       (persistence-lg-${real}   ${libraries})
set_target_properties       (persistence-lg-${real}   PROPERTIES COMPILE_FLAGS -DREEBER_REAL=${real})

add_executable              (persistent-integral-lg-${real}   persistent-integral-lg.cpp ${DEBUG_SOURCES} ${SOURCES})
target_link_libraries       (persistent-integral-lg-${real}   ${libraries})
if                          (persistent-integral-trace-vertices)
set_target_properties       (persistent-integral-lg-${real}   PROPERTIES COMPILE_FLAGS "-DREEBER_REAL=${real} -DREEBER_PERSISTENT_INTEGRAL_TRACE_VTCS")
else                        (persistent-integral-trace-vertices)
set_target_properties       (persistent-integral-lg-${real}   PROPERTIES COMPILE_FLAGS -DREEBER_REAL=${real})
endif                       (persistent-integral-trace-vertices)

add_executable              (nested-integral-lg-${real}   nested-integral-lg.cpp ${DEBUG_SOURCES} ${SOURCES})
target_link_libraries       (nested-integral-lg-${real}   ${libraries})
if                          (persistent-integral-trace-vertices)
set_target_properties       (nested-integral-lg-${real}   PROPERTIES COMPILE_FLAGS "-DREEBER_REAL=${real} -DREEBER_PERSISTENT_INTEGRAL_TRACE_VTCS")
else                        (persistent-integral-trace-vertices)
set_target_properties       (nested-integral-lg-${real}   PROPERTIES COMPILE_FLAGS -DREEBER_REAL=${real})
endif                       (persistent-integral-trace-vertices)

add_executable              (print-trees-${real}   print-trees.cpp ${DEBUG_SOURCES} ${SOURCES})
target_link_libraries       (print-trees-${real}   ${libraries})
set_target_properties       (print-trees-${real}   PROPERTIES COMPILE_FLAGS "-DREEBER_REAL=${real} -std=c++11")

endforeach                  (real)

if                          (boxlib)

string                      (TOLOWER ${CCSE_PRECISION} CCSE_PRECISION_LC)

target_link_libraries       (mt-local-global-${CCSE_PRECISION_LC}   ${libraries} ${CCSE_LIBRARIES})
set_target_properties       (mt-local-global-${CCSE_PRECISION_LC}   PROPERTIES COMPILE_FLAGS "-DREEBER_USE_BOXLIB_READER")

target_link_libraries       (mt-lg-ghosts-${CCSE_PRECISION_LC}      ${libraries} ${CCSE_LIBRARIES})
set_target_properties       (mt-lg-ghosts-${CCSE_PRECISION_LC}      PROPERTIES COMPILE_FLAGS "-DREEBER_USE_BOXLIB_READER")

target_link_libraries       (persistent-integral-lg-${CCSE_PRECISION_LC}  ${libraries} ${CCSE_LIBRARIES})

if                          (persistent-integral-trace-vertices)
set_target_properties       (persistent-integral-lg-${CCSE_PRECISION_LC}  PROPERTIES COMPILE_FLAGS "-DREEBER_USE_BOXLIB_READER -DREEBER_PERSISTENT_INTEGRAL_TRACE_VTCS")
else                        (persistent-integral-trace-vertices)
set_target_properties       (persistent-integral-lg-${CCSE_PRECISION_LC}  PROPERTIES COMPILE_FLAGS "-DREEBER_USE_BOXLIB_READER")
endif                       (persistent-integral-trace-vertices)
endif                       (boxlib)

add_subdirectory            (tests)
