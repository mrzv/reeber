# src files
foreach(real   float )

# AMR
find_package(AMReX REQUIRED HINTS ${AMREX_INSTALL_DIR})

set (libraries ${libraries}
     ${CMAKE_THREAD_LIBS_INIT}
     AMReX::amrex
     )

# include files
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
                    ${CMAKE_CURRENT_SOURCE_DIR}/../amr-merge-tree/include)

add_executable(amr_connected_components_${real}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/amr-connected-components.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../amr-merge-tree/src/amr-plot-reader.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../amr-merge-tree/src/hdf5-reader.cpp)

target_compile_definitions(amr_connected_components_${real} PRIVATE REEBER_REAL=${real} AMREX_SPACEDIM=3 BL_USE_MPI=ON ENABLE_DP=OFF BL_USE_FLOAT=1)

# target_link_libraries(amr_merge_tree PUBLIC ${libraries})
target_link_libraries(amr_connected_components_${real} PUBLIC ${libraries} HighFive)

endforeach()
