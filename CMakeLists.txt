cmake_minimum_required(VERSION 3.10)
project(adaptnotch)

set(CMAKE_CXX_STANDARD 14)
if(NOT CMAKE_BUILD_TYPE)
    # Options: Debug, Release, MinSizeRel, RelWithDebInfo
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

find_package(Eigen3 REQUIRED)
message(STATUS "Eigen Version: " ${EIGEN3_VERSION_STRING})

add_library(adaptnotch SHARED src/adaptnotch.cpp)
target_include_directories(adaptnotch PUBLIC include ${EIGEN3_INCLUDE_DIRS})

add_executable(main src/main.cpp)
target_include_directories(main PUBLIC ${EIGEN3_INCLUDE_DIRS})
target_link_libraries(main adaptnotch pthread)

add_executable(fft src/test.cpp)
target_include_directories(fft PUBLIC ${EIGEN3_INCLUDE_DIRS})
target_link_libraries(fft)