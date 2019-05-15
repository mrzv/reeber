cmake_minimum_required(VERSION 3.2)
project(amr_connected_components)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_definitions(-std=c++14)

# Default to Release
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif (NOT CMAKE_BUILD_TYPE)


set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g  -D_GLIBCXX_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 ")
#set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -O2 -g2 -debug full ")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -O2 -g -DNDEBUG -pg ")

# MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
set(libraries ${libraries} ${MPI_LIBRARIES})

# DIY
include_directories(${DIY_PATH}/include)

# threads
find_package(Threads)
set(libraries ${libraries} ${CMAKE_THREAD_LIBS_INIT})



# src files
foreach(real   float )

# AMR
if (DEFINED ENV{NERSC_HOST})
    set(AMREX_HOME "/global/homes/g/greynarn/code/amrex")
    message(STATUS "Looking for AMR on NERSC")
    if ("$ENV{NERSC_HOST}" STREQUAL "edison")
        message(STATUS "Linking with Edison version of AMReX")
        find_package(AMReX 18 REQUIRED HINTS "/global/homes/g/greynarn/code/amrex/install_edison_gnu_${real}")
    else()
        #message(STATUS "Linking with Cori (non-Edison) version of AMReX, intel compiled")
        #find_package(AMReX 18 REQUIRED HINTS "/global/homes/g/greynarn/code/amrex/install_cori_knl_intel_${real}")
        message(STATUS "Linking with Cori (non-Edison) version of AMReX, GNU compiled")
        find_package(AMReX 18 REQUIRED HINTS "/global/homes/g/greynarn/code/amrex/install_cori_knl_gnu_${real}")
    endif()
else()
    set(AMREX_HOME "/home/narn/code/berkeley-2018/amrex")
    find_package(AMReX 18 REQUIRED HINTS "/home/narn/code/berkeley-2018/amrex_install_dir_${real}")
    message(STATUS "Looking for AMR on Arnur's computer for ${real}")
endif()

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
set_target_properties(amr_connected_components_${real}   PROPERTIES COMPILE_DEFINITIONS "REEBER_REAL=${real}")

# target_link_libraries(amr_merge_tree PUBLIC ${libraries})
target_link_libraries(amr_connected_components_${real} PUBLIC ${libraries})

if (MPI_COMPILE_FLAGS)
    set_target_properties(amr_connected_components_${real} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif ()

if (MPI_LINK_FLAGS)
    set_target_properties(amr_connected_components_${real} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif ()

endforeach()
