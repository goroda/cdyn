cmake_minimum_required (VERSION 2.6)
project (CDYN)

#create my librray
set(C3_VERSION_MAJOR 1 CACHE STRING "major version" FORCE)
set(C3_VERSION_MINOR 0 CACHE STRING "minor version" FORCE)
set(C3_VERSION "${C3_VERSION_MAJOR}.${C3_VERSION_MINOR}" CACHE STRING "version" FORCE)

option(C3_DEBUG_MODE "whether or not debug mode is activated" OFF)

set(include_dest "include/cdyn")
set(main_lib_dest "lib/cdyn")
set(lib_dest "${main_lib_dest}/${CMAKE_BUILD_TYPE}")

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

    include_directories(/usr/include)
    set(blas_lapack "-lblas -llapack")

    #find_package(SWIG REQUIRED)
    #include(${SWIG_USE_FILE})
endif()

set(CMAKE_CXX_COMPILER "gcc")

#set(static_libs 1)
set(static_libs 0)

include_directories(${BLAS_INCLUDE_DIR}, ${LAPACK_INCLUDE_DIR})

set(CMAKE_C_FLAGS "${compiler_flags}")

#c3 library include files
#include_directories(../c3/include)

#include_directories(include)
include_directories(src)


add_subdirectory(src)
add_subdirectory(test)
#add_subdirectory(examples)


install(FILES cmake/cdyn-config.cmake DESTINATION ${main_lib_dest})
install(EXPORT cdyn DESTINATION ${lib_dest})
