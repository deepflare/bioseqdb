cmake_minimum_required(VERSION 3.14)

project(bioseqdb LANGUAGES CXX)

set (BIOSEQDB_ROOT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set (BIOSEQDB_DEPENDENCIES_PATH ${BIOSEQDB_ROOT_PATH}/dependencies) 
set (BIOSEQDB_PACKAGES_PATH ${BIOSEQDB_ROOT_PATH}/packages)
set (BIOSEQDB_COMMON_PATH ${BIOSEQDB_PACKAGES_PATH}/common)

add_subdirectory(dependencies/rpclib)

add_subdirectory(packages/bioseqdb_pg)
add_subdirectory(packages/mmseq2)