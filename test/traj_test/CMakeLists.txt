cmake_minimum_required (VERSION 2.6)

add_executable(traj_test traj_test.c)
TARGET_LINK_LIBRARIES(traj_test 
    -g 
    cdyn
    cutest 
    -lm)

