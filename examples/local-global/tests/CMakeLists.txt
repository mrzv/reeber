add_custom_target           (local-global-tests ALL
                             ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
                            )

add_test                    (dens40-pi          dens40-pi.sh)

include                     (ExternalData)
set                         (ExternalData_URL_TEMPLATES "http://data.mrzv.org/%(algo)/%(hash)")

if                          (slow-tests)
    ExternalData_Add_Test   (LyA-density-plt00242-128
                             NAME    LyAReeberTest-pi
                             COMMAND LyAReeberTest-pi.sh DATA{LyAReeberTest_plt00242_128x128x128_density.npy})
    ExternalData_Add_Target (LyA-density-plt00242-128)
endif                       (slow-tests)
