cmake_minimum_required(VERSION 3.0)
project(adaptnotch VERSION 0.1)

set(CMAKE_CXX_STANDARD 14)
if(NOT CMAKE_BUILD_TYPE)
    # Options: Debug, Release, MinSizeRel, RelWithDebInfo
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

find_package(Eigen3 REQUIRED)
message(STATUS "Eigen Version: " ${EIGEN3_VERSION_STRING})

add_library(adaptnotch SHARED src/adaptnotch.cpp)
target_include_directories(adaptnotch PUBLIC
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  ${EIGEN3_INCLUDE_DIRS})

option(BUILD_EXAPPS "Build adaptnotch example" ON)

if(BUILD_EXAPPS)
  # Setup cmake paths
  set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

  # Download plot library
  set(PLOT_DIR ${CMAKE_CURRENT_BINARY_DIR}/plot-download)
  set(BUILD_EXAMPLES OFF CACHE INTERNAL "") # don't build plot examples
  configure_file("${CMAKE_MODULE_PATH}/plot.cmake.in" "${PLOT_DIR}/CMakeLists.txt" IMMEDIATE @ONLY)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" . WORKING_DIRECTORY ${PLOT_DIR} )
  execute_process(COMMAND ${CMAKE_COMMAND} --build . WORKING_DIRECTORY ${PLOT_DIR})
  add_subdirectory(${PLOT_DIR}/src ${PLOT_DIR}/build)

  add_executable(main src/main.cpp)
  target_include_directories(main PUBLIC ${EIGEN3_INCLUDE_DIRS})
  target_link_libraries(main adaptnotch pthread plot)
endif()


include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/adaptnotch)

install(TARGETS adaptnotch
    EXPORT adaptnotch-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT adaptnotch-targets
  FILE adaptnotch-targets.cmake
  DESTINATION ${INSTALL_CONFIGDIR})

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/adaptnotch-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/adaptnotch-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/adaptnotch-config.cmake
  INSTALL_DESTINATION ${INSTALL_CONFIGDIR})

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/adaptnotch-config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/adaptnotch-config-version.cmake
  DESTINATION ${INSTALL_CONFIGDIR})