cmake_minimum_required(VERSION 3.14)

project(bioseqdb LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

set(PostgreSQL_ADDITIONAL_VERSIONS "14")
set (POSTGRES_VERSION 14)
set(CMAKE_VERBOSE_MAKEFILE ON)

add_definitions(-DSIMDE_ENABLE_NATIVE_ALIASES)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/dependencies/simde")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/dependencies")


find_package(PostgreSQL ${POSTGRES_VERSION} EXACT REQUIRED)
find_program(PG_CONFIG pg_config)
execute_process(COMMAND ${PG_CONFIG} --pkglibdir OUTPUT_VARIABLE PG_CONFIG_PKGLIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${PG_CONFIG} --sharedir OUTPUT_VARIABLE PG_CONFIG_SHAREDIR OUTPUT_STRIP_TRAILING_WHITESPACE)

find_package(BZip2 REQUIRED)
include(ExternalProject)

set(EXTERNAL_C_FLAGS "-g -Wall -O2 -fPIC -march=native -I${CMAKE_CURRENT_SOURCE_DIR}/dependencies/simde -DSIMDE_ENABLE_NATIVE_ALIASES")

ExternalProject_Add(bwamem
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/bwa"
  BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/bwamem"
  BUILD_COMMAND make CFLAGS=${EXTERNAL_C_FLAGS} CXXFLAGS=${EXTERNAL_C_FLAGS} libbwa.a
  BUILD_IN_SOURCE 0
  CONFIGURE_COMMAND cp -R ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/bwa/. .
  UPDATE_COMMAND "" # Skip annoying updates for every build
  INSTALL_COMMAND ""
)

ExternalProject_Get_Property(bwamem install_dir)

add_library(libbwa SHARED IMPORTED)
set_target_properties(libbwa PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/bwamem/libbwa.a)
add_dependencies(libbwa bwamem)

add_library(bioseqdb_pg SHARED
        bioseqdb_pg/bwa.cpp
        bioseqdb_pg/extension.cpp
        bioseqdb_pg/sequence.cpp
        )
add_executable(bioseqdb_import
        bioseqdb_import/main.cpp
        )

add_dependencies(bioseqdb_pg libbwa)

target_include_directories(bioseqdb_pg PRIVATE ${PostgreSQL_TYPE_INCLUDE_DIR})
target_link_libraries(bioseqdb_pg PRIVATE ${PostgreSQL_LIBRARIES} ${BZIP2_LIBRARIES} libbwa)


target_include_directories(bioseqdb_import PRIVATE ${PostgreSQL_INCLUDE_DIRS})
target_link_libraries(bioseqdb_import PRIVATE ${PostgreSQL_LIBRARIES})

install(TARGETS bioseqdb_pg DESTINATION ${PG_CONFIG_PKGLIBDIR})
install(FILES bioseqdb_pg/bioseqdb.control DESTINATION ${PG_CONFIG_SHAREDIR}/extension)
install(FILES bioseqdb_pg/bioseqdb--0.0.0.sql DESTINATION ${PG_CONFIG_SHAREDIR}/extension)
