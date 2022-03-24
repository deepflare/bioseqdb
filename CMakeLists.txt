cmake_minimum_required(VERSION 3.14)

project(bioseqdb LANGUAGES CXX)

find_package(PostgreSQL REQUIRED COMPONENTS Server)

add_library(bioseqdb_pg SHARED
        bioseqdb_pg/extension.cpp
        )

target_include_directories(bioseqdb_pg PRIVATE ${PostgreSQL_TYPE_INCLUDE_DIR})
target_link_libraries(bioseqdb_pg PRIVATE ${PostgreSQL_LIBRARIES})
