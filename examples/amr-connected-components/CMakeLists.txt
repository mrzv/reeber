# src files
foreach(real   float )

# AMR
find_package(AMReX 18 REQUIRED HINTS ${AMREX_INSTALL_DIR})

set (libraries ${libraries}
        ${CMAKE_THREAD_LIBS_INIT}
        AMReX::amrex
        gfortran
        )

# include files
file(GLOB AMT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp)
file(GLOB REEBER_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/../../include/reeber/*.h ${CMAKE_CURRENT_SOURCE_DIR}/../../include/reeber/*.hpp)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
                    ${CMAKE_CURRENT_SOURCE_DIR}/../amr-merge-tree/include
                    ${AMREX_HOME}/Src/Base
                    ${AMREX_HOME}/Src/Extern/amrdata)

add_executable(amr_connected_components_${real}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/amr-connected-components.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../amr-merge-tree/src/amr-plot-reader.cpp
        ${AMT_HEADERS}
        ${REEBER_HEADERS})

target_compile_definitions(amr_connected_components_${real} PRIVATE REEBER_REAL=${real} AMREX_SPACEDIM=3 BL_USE_MPI=ON ENABLE_DP=OFF BL_USE_FLOAT=1)

# target_link_libraries(amr_merge_tree PUBLIC ${libraries})
target_link_libraries(amr_connected_components_${real} PUBLIC ${libraries})

if (MPI_COMPILE_FLAGS)
    set_target_properties(amr_connected_components_${real} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif ()

if (MPI_LINK_FLAGS)
    set_target_properties(amr_connected_components_${real} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif ()

endforeach()
