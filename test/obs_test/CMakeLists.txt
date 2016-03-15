cmake_minimum_required (VERSION 2.6)

add_executable(obs_test obs_test.c)
TARGET_LINK_LIBRARIES(obs_test 
    -g 
    cdyn
    cutest 
    -lm)

