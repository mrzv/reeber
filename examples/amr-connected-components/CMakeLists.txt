# src files
foreach(real   float )

set (libraries ${libraries}
     ${CMAKE_THREAD_LIBS_INIT}
     )

# include files
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
                    ${CMAKE_CURRENT_SOURCE_DIR}/../amr-merge-tree/include)

set(amr_cc_sources  ${CMAKE_CURRENT_SOURCE_DIR}/src/amr-connected-components.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/../amr-merge-tree/src/hdf5-reader.cpp)


# AMR
if (amrex)
    find_package(AMReX REQUIRED HINTS ${AMREX_INSTALL_DIR})
    set (libraries ${libraries} AMReX::amrex)
    set(amr_cc_sources  ${amr_cc_sources} ${CMAKE_CURRENT_SOURCE_DIR}/../amr-merge-tree/src/amr-plot-reader.cpp)
endif()

add_executable(amr_connected_components_${real} ${amr_cc_sources})

target_compile_definitions(amr_connected_components_${real} PRIVATE REEBER_REAL=${real} AMREX_SPACEDIM=3 BL_USE_MPI=ON ENABLE_DP=OFF BL_USE_FLOAT=1)

# target_link_libraries(amr_merge_tree PUBLIC ${libraries})
target_link_libraries(amr_connected_components_${real} PUBLIC ${libraries} HighFive)

endforeach()
