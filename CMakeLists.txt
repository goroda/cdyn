cmake_minimum_required (VERSION 2.6)
project (CDYN)

#create my librray
set(CDYN_VERSION_MAJOR 1 CACHE STRING "major version" FORCE)
set(CDYN_VERSION_MINOR 0 CACHE STRING "minor version" FORCE)
set(CDYN_VERSION "${CDYN_VERSION_MAJOR}.${CDYN_VERSION_MINOR}" CACHE STRING "version" FORCE)

option(CDYN_DEBUG_MODE "whether or not debug mode is activated" OFF)

if (APPLE)
    option(USE_CLANG "build application with clang" ON)
    include_directories(/opt/local/include)
    find_library(macblas
        NAMES ACCELERATE
        PATHS ${CMAKE_OSX_SYSROOT}/System/Library
        PATH_SUFFIXES Frameworks
        NO_DEFAULT_PATH)
        
    set(blas_lapack ${macblas}/Accelerate)
    set(compiler_flags "-std=c99 -O3 -Wall -Weverything -Wno-shorten-64-to-32 -Wno-padded -Wno-missing-prototypes -Wno-format-nonliteral -Wno-disabled-macro-expansion")
    set(CMAKE_C_COMPILER "clang")
else ()
    set(compiler_flags "-g -Wall -Wextra -O3 -std=c99")
    set(CMAKE_C_COMPILER "gcc")

    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    include_directories(${BLAS_INCLUDE_DIR}, ${LAPACK_INCLUDE_DIR})
    set(blas_lapack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    
    include_directories(/usr/include)
#   set(blas_lapack "-lblas -llapack")
endif()

set(CMAKE_CXX_COMPILER "gcc")

set(static_libs 0)

include_directories(${BLAS_INCLUDE_DIR}, ${LAPACK_INCLUDE_DIR})

set(CMAKE_C_FLAGS "${compiler_flags}")

include_directories(src)
add_subdirectory(src)
add_subdirectory(test)

