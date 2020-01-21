cmake_minimum_required(VERSION 2.8)
project(tv_min_test)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

add_subdirectory(googletest/googletest)

add_executable(tv_min_test tv_min_test.cpp )

target_link_libraries(tv_min_test 
						PRIVATE gtest)