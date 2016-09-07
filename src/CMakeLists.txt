cmake_minimum_required (VERSION 2.6)

file(GLOB SRC_FILES dynamics.c observations.c integrate.c simulate.c)

set(header integrate.h
           observations.h
           simulate.h)
           
set(src observations.c
        integrate.c
        simulate.c)


# if ( ${static_libs} )
#     add_library(cdyn ${SRC_FILES})
# else ()
#     add_library(cdyn SHARED ${SRC_FILES})
#     target_link_libraries(cdyn
#         -fPIC
#         -lm
#         ${blas_lapack})
# endif() 

add_library(cdyn SHARED ${header} ${src})
target_link_libraries(cdyn
  ${blas_lapack}
  -lm
  )
target_include_directories(cdyn PUBLIC 
                     $<BUILD_INTERFACE:${CDYN_SOURCE_DIR}/include>
                     $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
                     $<INSTALL_INTERFACE:${include_dest}>
                     $<INSTALL_INTERFACE:${lib_dest}>)

install(TARGETS cdyn EXPORT cdyn DESTINATION "${lib_dest}")
install(FILES ${header} DESTINATION "${include_dest}")

