cmake_minimum_required (VERSION 2.6)

add_executable(integrate_test integrate_test.c)
TARGET_LINK_LIBRARIES(integrate_test 
    -g 
    cdyn
    cutest 
    -lm)

