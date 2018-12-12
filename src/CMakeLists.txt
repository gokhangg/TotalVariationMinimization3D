cmake_minimum_required(VERSION 2.8)
project(tv_min)

add_executable(tv_min tv_min.cpp )

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
target_link_libraries(tv_min ${ITK_LIBRARIES})
