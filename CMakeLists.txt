cmake_minimum_required(VERSION 3.14)

project(bioseqdb LANGUAGES CXX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic")

add_custom_target(
    my_custom_target_that_always_runs ALL
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/__header.h
)

add_custom_command(
    OUTPUT
        ${CMAKE_CURRENT_BINARY_DIR}/__header.h  # fake! ensure we run!
    # this command must generate: ${CMAKE_CURRENT_BINARY_DIR}/header.h
    COMMAND poetry install
)


set (BIOSEQDB_ROOT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set (BIOSEQDB_DEPENDENCIES_PATH ${BIOSEQDB_ROOT_PATH}/dependencies) 
set (BIOSEQDB_PACKAGES_PATH ${BIOSEQDB_ROOT_PATH}/packages)
set (BIOSEQDB_COMMON_PATH ${BIOSEQDB_PACKAGES_PATH}/common)

set(POSTGRES_EXTENSION_NAME bioseqdb)
set(POSTGRES_EXTENSION_VERSION 0.0.0)

add_subdirectory(packages/bioseqdb_pg)
add_subdirectory(dependencies/rpclib)
add_subdirectory(packages/mmseq2)