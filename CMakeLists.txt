project                     (Reeber Fortran CXX C)
cmake_minimum_required      (VERSION 3.15)

option                      (amrex              "Build Reeber with AMReX support"               OFF)
option                      (trace              "Enable trace log"                              OFF)
option                      (profile            "Build Reeber with profiling"                   ON)
option                      (counters           "Build Reeber with counters"                    OFF)
option                      (slow-tests         "Enable slow tests"                             ON)
option                      (use-tbb            "Thread using TBB"                              OFF)

add_definitions             (-Wall -fPIC)

# Default to Release
if                          (NOT CMAKE_BUILD_TYPE)
    set                     (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
    set_property            (CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif                       (NOT CMAKE_BUILD_TYPE)

# Debugging
if                          (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    add_definitions         (-DDEBUG)
endif                       ()

# Find everything that's always required
find_package                (Boost)

find_package                (MPI REQUIRED)

find_package                (Threads)

# TBB
if                          (use-tbb)
find_path                   (TBB_INCLUDE_DIRS tbb/tbb.h)
find_library                (TBB_LIBRARY NAMES tbb)
find_library                (TBB_MALLOC_LIBRARY NAMES tbbmalloc)

   if                       (TBB_INCLUDE_DIRS AND TBB_LIBRARY AND TBB_MALLOC_LIBRARY)
      add_definitions       (-DREEBER_USE_TBB)
      include_directories   (${TBB_INCLUDE_DIRS})
      set                   (CMAKE_CXX_FLAGS "-mcx16 ${CMAKE_CXX_FLAGS}")
   else                     (TBB_FOUND)
      message               ("TBB not found; disabling")
   endif                    (TBB_INCLUDE_DIRS AND TBB_LIBRARY AND TBB_MALLOC_LIBRARY)
endif                       (use-tbb)

set                         (CMAKE_CXX_STANDARD 14)

if                          (profile)
    add_definitions         (-DPROFILE)
endif                       (profile)

# Trace logging (won't do much without debugging being on)
if                          (trace)
    add_definitions         (-DTRACE)
endif                       (trace)

# Counters
if                          (counters)
    add_definitions         (-DCOUNTERS)
endif                       (counters)

# DIY
find_path                   (DIY_INCLUDE_DIR        diy/master.hpp)

# HighFive
add_subdirectory            (ext/HighFive)

# Set includes
include_directories         (${CMAKE_CURRENT_SOURCE_DIR}/include
                             SYSTEM ${DIY_INCLUDE_DIR}
                             SYSTEM ${MPI_INCLUDE_PATH}
                             SYSTEM ${Boost_INCLUDE_DIR})

# Set libraries
set                         (libraries
                             ${CMAKE_THREAD_LIBS_INIT}
                             ${TBB_LIBRARY}
                             ${TBB_MALLOC_LIBRARY}
                             ${MPI_C_LIBRARIES}
                             ${MPI_CXX_LIBRARIES})

# backwards.hpp
if                          (debug)
    find_library            (LIBDW_LIBRARY NAMES dw)
    if                      (LIBDW_LIBRARY)
        set                 (DEBUG_SOURCES ${CMAKE_SOURCE_DIR}/src/backward.cpp)
        add_definitions     (-DBACKWARD_HAS_DW=1)
        set                 (libraries ${libraries} ${LIBDW_LIBRARY})
    else                    (LIBDW_LIBRARY)
        message             (STATUS "LibDW not found; backward.cpp won't be used")
    endif                   (LIBDW_LIBRARY)
endif                       (debug)

enable_testing              ()
add_subdirectory            (examples)
