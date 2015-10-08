cmake_minimum_required(VERSION 2.8.4)

#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

project(compsph)

find_package( OpenCV REQUIRED )
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES 
	compsph.cpp 
	../yuv/yuv_helper.cpp 
	../map_utils/map_utils.cpp
	../panomapper/panomapper.cpp
	)

include_directories(../yuv ../map_utils ../panomapper)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${compsph_SOURCE_DIR})

add_executable(compsph ${SOURCE_FILES})

target_link_libraries(compsph ${OpenCV_LIBS} )
