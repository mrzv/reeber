add_custom_target           (local-global-tests ALL
                             ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
                            )

add_test                    (dens40-pi          dens40-pi.sh)
