cmake_minimum_required (VERSION 2.6)

include_directories(CuTest)
add_subdirectory(CuTest)

add_subdirectory(dyn_test)
add_subdirectory(obs_test)
add_subdirectory(int_test)
add_subdirectory(traj_test)


