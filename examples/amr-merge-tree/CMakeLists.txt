#foreach(real   float double)
foreach(real   float)
# AMR
find_package(AMReX REQUIRED HINTS ${AMREX_INSTALL_DIR})

set (libraries ${libraries}
        ${CMAKE_THREAD_LIBS_INIT}
        AMReX::amrex
        gfortran            # ????
        )

# includes
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
        ${DIY_PATH}/include
        ${AMREX_HOME}/Src/Base
        ${AMREX_HOME}/Src/Extern/amrdata
        #${AMREX_HOME}/Src/AmrCore
        )

add_executable(amr_merge_tree_simple_${real} ${CMAKE_CURRENT_SOURCE_DIR}/src/amr-merge-tree.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/amr-plot-reader.cpp)
if ("${real}" STREQUAL "float")
    target_compile_definitions(amr_merge_tree_simple_${real} PRIVATE REEBER_REAL=${real} AMREX_SPACEDIM=3 BL_USE_MPI=ON ENABLE_DP=OFF BL_USE_FLOAT=1)
else()
    target_compile_definitions(amr_merge_tree_simple_${real} PRIVATE REEBER_REAL=${real} AMREX_SPACEDIM=3 BL_USE_MPI=ON ENABLE_DP=ON BL_USE_FLOAT=0)
#    ket_target_properties(amr_merge_tree_simple_${real}   PROPERTIES COMPILE_DEFINITIONS "REEBER_REAL=${real};ENABLE_DP=ON;BL_USE_DOUBLE=1")
endif()


add_executable(amr_merge_tree_test_${real} ${CMAKE_CURRENT_SOURCE_DIR}/tests/tests_main.cpp ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_amr_merge_tree.cpp)
set_target_properties(amr_merge_tree_test_${real} PROPERTIES COMPILE_DEFINITIONS "REEBER_REAL=${real}")

add_executable(write_refined_amr_${real} ${CMAKE_CURRENT_SOURCE_DIR}/src/write-refined-amr.cpp)
set_target_properties(write_refined_amr_${real}  PROPERTIES COMPILE_DEFINITIONS "REEBER_REAL=${real}")


target_link_libraries(amr_merge_tree_simple_${real} PUBLIC ${libraries})
target_link_libraries(amr_merge_tree_test_${real} PUBLIC ${libraries})
target_link_libraries(write_refined_amr_${real} PUBLIC ${libraries})

add_test(amr-merge-tree-test-${real} amr_merge_tree_test_${real})

endforeach()
