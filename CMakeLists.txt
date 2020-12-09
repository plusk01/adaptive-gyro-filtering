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

# Setup cmake paths
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Download plot library
set(PLOT_DIR ${CMAKE_CURRENT_BINARY_DIR}/plot-download)
set(BUILD_EXAMPLES OFF CACHE BOOL "") # don't build plot examples
configure_file("${CMAKE_MODULE_PATH}/plot.cmake.in" "${PLOT_DIR}/CMakeLists.txt" IMMEDIATE @ONLY)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" . WORKING_DIRECTORY ${PLOT_DIR} )
execute_process(COMMAND ${CMAKE_COMMAND} --build . WORKING_DIRECTORY ${PLOT_DIR})
add_subdirectory(${PLOT_DIR}/src ${PLOT_DIR}/build)


add_library(adaptnotch SHARED src/adaptnotch.cpp)
target_include_directories(adaptnotch PUBLIC include ${EIGEN3_INCLUDE_DIRS})


add_executable(main src/main.cpp)
target_include_directories(main PUBLIC ${EIGEN3_INCLUDE_DIRS})
target_link_libraries(main adaptnotch pthread plot)
