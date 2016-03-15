cmake_minimum_required (VERSION 2.6)

file(GLOB SRC_FILES dynamics.c observations.c integrate.c simulate.c)

if ( ${static_libs} )
    add_library(cdyn ${SRC_FILES})
else ()
    add_library(cdyn SHARED ${SRC_FILES})
    target_link_libraries(cdyn
        -fPIC
        -lm
        ${blas_lapack})
endif() 

