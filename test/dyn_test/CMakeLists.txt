cmake_minimum_required (VERSION 2.6)

add_executable(dyn_test dyn_test.c)
TARGET_LINK_LIBRARIES(dyn_test 
    -g 
    cdyn
    cutest 
    -lm)

