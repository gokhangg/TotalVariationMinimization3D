cmake_minimum_required(VERSION 3.0)
project(TV_MIN)

#set(ITK_DIR XXX)
set(DATASET_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Dataset)
set(TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Test)
set(PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

option(BUILD_TEST "Select this if you want tests for the TVimage" OFF)
option(DOWNLOAD_TEST_DATA "Download a test image" OFF)

add_subdirectory(${PROJECT_SOURCE_DIR})

if (BUILD_TEST)
	add_subdirectory(${TEST_DIR})
endif()

if (DOWNLOAD_TEST_DATA)
	if (NOT EXISTS ${DATASET_DIR})
		execute_process(COMMAND git clone https://github.com/gokhangg/Datasets.git ${DATASET_DIR})
	endif()
endif()