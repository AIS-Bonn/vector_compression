
cmake_minimum_required(VERSION 2.8)
project(vector_compression)

# Detect if we are building under catkin
find_package(catkin QUIET)

if(catkin_FOUND)
	catkin_package(
		INCLUDE_DIRS include
		LIBRARIES vector_compression
	)
endif()

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

include_directories(include)

add_definitions("-std=c++11" "-Wall")

add_library(vector_compression SHARED
	src/vector_compression.cpp
)

if(NOT catkin_FOUND)
	install(
		TARGETS vector_compression
		LIBRARY DESTINATION lib
	)
	install(
		DIRECTORY include/vector_compression
		DESTINATION include
	)
endif()

add_executable(vector_compression_tests
	tests/test_main.cpp
	tests/test_vector_compression.cpp
)
target_link_libraries(vector_compression_tests
	vector_compression
)
