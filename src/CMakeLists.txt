cmake_minimum_required(VERSION 3.0)
project(TV_MIN_FILTER)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

set(CURRENT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(TVIMAGE_DIR ${CURRENT_DIR}/TV_Image)
set(COMMANDPARSER_DIR ${CMAKE_BINARY_DIR}/CommandlineParser)

if (NOT EXISTS ${COMMANDPARSER_DIR})
	execute_process(COMMAND git clone https://github.com/gokhangg/CommandlineParser.git ${COMMANDPARSER_DIR})
endif()

include_directories(${TVIMAGE_DIR} ${COMMANDPARSER_DIR})
set(HEADER_FILES tv_filter.h tv_filter.hxx ${TVIMAGE_DIR}/tv_image.h)
add_executable(TV_MIN_FILTER TV_MIN.cpp ${HEADER_FILES})
target_link_libraries(TV_MIN_FILTER ${ITK_LIBRARIES})
