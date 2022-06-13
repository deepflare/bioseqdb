project(HelloWorld)
include(ExternalProject)
cmake_minimum_required(VERSION 3.0)

add_definitions(-DSIMDE_ENABLE_NATIVE_ALIASES)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/simde")

ExternalProject_Add(bwa
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/bwa"
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/bwa
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make libbwa.a
  BUILD_IN_SOURCE 1
)

add_executable(hello_world main.cpp)