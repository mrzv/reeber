cmake_minimum_required(VERSION 3.2)
project(amr_merge_tree)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_definitions(-std=c++14)

# Default to Release
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif (NOT CMAKE_BUILD_TYPE)


set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb  -D_GLIBCXX_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 ")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -O2 -g -ggdb -pg")

# MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
set(libraries ${libraries} ${MPI_LIBRARIES})

# DIY
include_directories(${DIY_PATH}/include)

# threads
find_package(Threads)
set(libraries ${libraries} ${CMAKE_THREAD_LIBS_INIT})


# AMR
set(AMREX_HOME "/home/narn/code/berkeley-2018/amrex")
find_package(AMReX 18 REQUIRED HINTS /home/narn/berkeley-2018/amrex_install_dir_14/)
#
#
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

# src files
#file(GLOB AMT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)


foreach(real   float double)

add_executable(amr_merge_tree_simple_${real} ${CMAKE_CURRENT_SOURCE_DIR}/src/amr-merge-tree.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/amr-plot-reader.cpp ${AMT_HEADERS} ${REEBER_HEADERS})
set_target_properties(amr_merge_tree_simple_${real}   PROPERTIES COMPILE_DEFINITIONS "REEBER_REAL=${real}")

#add_executable(amr_merge_tree_send_components_${real} ${CMAKE_CURRENT_SOURCE_DIR}/src/amr-merge-tree.cpp ${AMT_HEADERS} ${REEBER_HEADERS})
#set_target_properties(amr_merge_tree_send_components_${real} PROPERTIES COMPILE_DEFINITIONS "REEBER_REAL=${real};AMR_MT_SEND_COMPONENTS")

add_executable(amr_merge_tree_test_${real} ${CMAKE_CURRENT_SOURCE_DIR}/tests/tests_main.cpp ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_amr_merge_tree.cpp ${AMT_HEADERS} ${REEBER_HEADERS})
set_target_properties(amr_merge_tree_test_${real} PROPERTIES COMPILE_DEFINITIONS "REEBER_REAL=${real}")

add_executable(write_refined_amr_${real} ${CMAKE_CURRENT_SOURCE_DIR}/src/write-refined-amr.cpp ${AMT_HEADERS} ${REEBER_HEADERS})
set_target_properties(write_refined_amr_${real}  PROPERTIES COMPILE_DEFINITIONS "REEBER_REAL=${real}")

#add_executable(my_tests ${CMAKE_CURRENT_SOURCE_DIR}/src/tests.cpp ${AMT_HEADERS} ${REEBER_HEADERS})

target_link_libraries(amr_merge_tree_simple_${real} PUBLIC ${libraries})
#target_link_libraries(amr_merge_tree_send_components_${real} PUBLIC ${libraries})
target_link_libraries(amr_merge_tree_test_${real} PUBLIC ${libraries})
target_link_libraries(write_refined_amr_${real} PUBLIC ${libraries})

if (MPI_COMPILE_FLAGS)
    set_target_properties(amr_merge_tree_simple_${real} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
#    set_target_properties(amr_merge_tree_send_components_${real} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    set_target_properties(write_refined_amr_${real} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
#    set_target_properties(my_tests PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif ()

if (MPI_LINK_FLAGS)
    set_target_properties(amr_merge_tree_simple_${real} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
#    set_target_properties(amr_merge_tree_send_components_${real} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
    set_target_properties(write_refined_amr_${real} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
#    set_target_properties(my_tests PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif ()

add_definitions(-std=c++14)
endforeach()
