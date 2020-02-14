cmake_minimum_required(VERSION 3.0)
project(TV_IMAGE_TEST)

include(ExternalProject)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(BUILD_SHARED ON CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)

set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/gtestExternal)

ExternalProject_Add(googletest
    GIT_REPOSITORY https://github.com/google/googletest
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION}
)

include_directories(${EXTERNAL_INSTALL_LOCATION}/include ${CMAKE_SOURCE_DIR}/src/TV_Image)
link_directories(${EXTERNAL_INSTALL_LOCATION}/lib)

if (MSVC)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
endif()

add_executable(TV_IMAGE_TEST tv_image_test.cpp)

add_dependencies(TV_IMAGE_TEST googletest)
target_link_libraries(TV_IMAGE_TEST debug gtestd)
target_link_libraries(TV_IMAGE_TEST optimized gtest)
