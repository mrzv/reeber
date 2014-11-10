project                     (Reeber)
cmake_minimum_required      (VERSION 2.8)

option                      (trace              "Enable trace log"                              OFF)
option                      (debug              "Build Reeber with debugging on"                OFF)
option                      (optimize           "Build Reeber with optimization"                ON)
option                      (profile            "Build Reeber with profiling"                   ON)
option                      (counters           "Build Reeber with counters"                    OFF)

add_definitions             (-Wall -fPIC)

# Find everything that's always required
find_package                (Boost)
find_package                (MPI REQUIRED)

# Debugging
if                          (debug)
    if                      (optimize)
            set             (cxx_flags                  ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
    else                    (optimize)
            set             (cxx_flags                  ${CMAKE_CXX_FLAGS_DEBUG})
    endif                   (optimize)
    add_definitions         (-DDEBUG)
else                        (debug)
    if                      (optimize)
            set             (cxx_flags                  ${CMAKE_CXX_FLAGS_RELEASE})
    else                    (optimize)
            set             (cxx_flags                  ${CMAKE_CXX_FLAGS})
    endif                   (optimize)
endif                       (debug)
add_definitions             (${cxx_flags})

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

# Set includes
include_directories         (${CMAKE_CURRENT_SOURCE_DIR}/include
                             ${DIY_INCLUDE_DIR}
                             SYSTEM ${MPI_INCLUDE_PATH}
                             SYSTEM ${Boost_INCLUDE_DIR})

# Set libraries
set                         (libraries
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

add_subdirectory            (examples)
