cmake_minimum_required(VERSION 3.14)

project(bioseqdb LANGUAGES CXX)

set (BIOSEQDB_ROOT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set (BIOSEQDB_DEPENDENCIES_PATH ${BIOSEQDB_ROOT_PATH}/dependencies) 

add_subdirectory(packages/bioseqdb_pg)