#foreach(real   float double)
foreach(real   float)
# AMR
find_package(AMReX 18 REQUIRED HINTS ${AMREX_INSTALL_DIR})

#if (DEFINED ENV{NERSC_HOST})
#    set(AMREX_HOME "/global/homes/g/greynarn/code/amrex")
#    message(STATUS "Looking for AMR on NERSC")
#    if ("$ENV{NERSC_HOST}" STREQUAL "edison")
#        message(STATUS "Linking with Edison version of AMReX")
#        find_package(AMReX 18 REQUIRED HINTS "/global/homes/g/greynarn/code/amrex/install_edison_gnu_${real}")
#    else()
#        #message(STATUS "Linking with Cori (non-Edison) version of AMReX, intel")
#        #find_package(AMReX 18 REQUIRED HINTS "/global/homes/g/greynarn/code/amrex/install_cori_knl_intel_${real}")
#        message(STATUS "Linking with Cori (non-Edison) version of AMReX, GNU")
#        find_package(AMReX 18 REQUIRED HINTS "/global/homes/g/greynarn/code/amrex/install_cori_knl_gnu_${real}")
#    endif()
#else()
#    set(AMREX_HOME "/home/narn/code/berkeley-2019/amrex")
#    find_package(AMReX 18 REQUIRED HINTS "/home/narn/code/berkeley-2019/amrex/install")
#    message(STATUS "Looking for AMR on Arnur's computer for ${real}")
#endif()

set (libraries ${libraries}
        ${CMAKE_THREAD_LIBS_INIT}
        AMReX::amrex
        gfortran
        )

# include files
file(GLOB AMT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp)
file(GLOB REEBER_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/../../include/reeber/*.h ${CMAKE_CURRENT_SOURCE_DIR}/../../include/reeber/*.hpp)
# includes
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
        ${DIY_PATH}/include
        ${AMREX_HOME}/Src/Base
        ${AMREX_HOME}/Src/Extern/amrdata
        #${AMREX_HOME}/Src/AmrCore
        )

add_executable(amr_merge_tree_simple_${real} ${CMAKE_CURRENT_SOURCE_DIR}/src/amr-merge-tree.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/amr-plot-reader.cpp ${AMT_HEADERS} ${REEBER_HEADERS})
if ("${real}" STREQUAL "float")
    target_compile_definitions(amr_merge_tree_simple_${real} PRIVATE REEBER_REAL=${real} AMREX_SPACEDIM=3 BL_USE_MPI=ON ENABLE_DP=OFF BL_USE_FLOAT=1)
else()
    target_compile_definitions(amr_merge_tree_simple_${real} PRIVATE REEBER_REAL=${real} AMREX_SPACEDIM=3 BL_USE_MPI=ON ENABLE_DP=ON BL_USE_FLOAT=0)
#    ket_target_properties(amr_merge_tree_simple_${real}   PROPERTIES COMPILE_DEFINITIONS "REEBER_REAL=${real};ENABLE_DP=ON;BL_USE_DOUBLE=1")
endif()


add_executable(amr_merge_tree_test_${real} ${CMAKE_CURRENT_SOURCE_DIR}/tests/tests_main.cpp ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_amr_merge_tree.cpp ${AMT_HEADERS} ${REEBER_HEADERS})
set_target_properties(amr_merge_tree_test_${real} PROPERTIES COMPILE_DEFINITIONS "REEBER_REAL=${real}")

add_executable(write_refined_amr_${real} ${CMAKE_CURRENT_SOURCE_DIR}/src/write-refined-amr.cpp ${AMT_HEADERS} ${REEBER_HEADERS})
set_target_properties(write_refined_amr_${real}  PROPERTIES COMPILE_DEFINITIONS "REEBER_REAL=${real}")


target_link_libraries(amr_merge_tree_simple_${real} PUBLIC ${libraries})
target_link_libraries(amr_merge_tree_test_${real} PUBLIC ${libraries})
target_link_libraries(write_refined_amr_${real} PUBLIC ${libraries})

if (MPI_COMPILE_FLAGS)
    set_target_properties(amr_merge_tree_simple_${real} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    set_target_properties(write_refined_amr_${real} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif ()

if (MPI_LINK_FLAGS)
    set_target_properties(amr_merge_tree_simple_${real} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
    set_target_properties(write_refined_amr_${real} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif ()

endforeach()