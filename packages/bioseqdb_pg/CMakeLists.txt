cmake_minimum_required(VERSION 3.17)

project(bioseqdb LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic")

if (NOT DEFINED MMSEQ_HOSTNAME)
    set(MMSEQ_HOSTNAME "localhost")
endif()
if (NOT DEFINED MMSEQ_PORT)
    set(MMSEQ_PORT "8080")
endif()
if (NOT DEFINED MMSEQ_TIMEOUT)
    set(MMSEQ_TIMEOUT "30000")
endif()

set(PostgreSQL_ADDITIONAL_VERSIONS "14")
set (POSTGRES_VERSION 14)
set(CMAKE_VERBOSE_MAKEFILE ON)

add_definitions(-DSIMDE_ENABLE_NATIVE_ALIASES)
include_directories("${BIOSEQDB_DEPENDENCIES_PATH}/simde")
include_directories("${BIOSEQDB_DEPENDENCIES_PATH}")
include_directories("${BIOSEQDB_PACKAGES_PATH}")
include_directories(${BIOSEQDB_DEPENDENCIES_PATH}/rpclib/include)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/extension")

file(GLOB SRC_TEMPLATE_FILES ${PROJECT_SOURCE_DIR}/extension/*.jinja)
file(GLOB SRC_EXTENSION_FILES ${PROJECT_SOURCE_DIR}/extension/*/*)


function(print_directory_variables dir skip_regex)
    message(STATUS "Call print_directory_variables")  
    # Dump variables:
    get_property(_variableNames DIRECTORY ${dir} PROPERTY VARIABLES)
    list (SORT _variableNames)
    foreach (_variableName ${_variableNames})
          unset(MATCHED)
          string(REGEX MATCH CMAKE_CXX MATCHED ${_variableName})
          if (MATCHED)
              continue()
          endif()
        get_directory_property(_variableValue DIRECTORY ${dir} DEFINITION ${_variableName})
        set(VARS_LIST ${VARS_LIST} "${_variableName}='${_variableValue}'")
    endforeach()

    message(STATUS "VARS_LIST ${dir}: ${VARS_LIST}")

    find_program(CHMOD_PROG chmod)
    execute_process(COMMAND ${CHMOD_PROG} 771 ${CMAKE_CURRENT_SOURCE_DIR}/generate_sources.sh)

    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/gen/config.h ${CMAKE_CURRENT_SOURCE_DIR}/gen/extension.cpp ${CMAKE_CURRENT_SOURCE_DIR}/gen/${POSTGRES_EXTENSION_NAME}--${POSTGRES_EXTENSION_VERSION}.sql ${CMAKE_CURRENT_SOURCE_DIR}/gen/${POSTGRES_EXTENSION_NAME}.control
      COMMAND ${CMAKE_COMMAND} -E env ${VARS_LIST} ${CMAKE_CURRENT_SOURCE_DIR}/generate_sources.sh
      DEPENDS ${SRC_TEMPLATE_FILES} ${SRC_EXTENSION_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/generate_sources.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endfunction(print_directory_variables)

print_directory_variables(. CMAKE_CXX)

find_package(PostgreSQL ${POSTGRES_VERSION} EXACT REQUIRED)
find_program(PG_CONFIG pg_config)
execute_process(COMMAND ${PG_CONFIG} --pkglibdir OUTPUT_VARIABLE PG_CONFIG_PKGLIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${PG_CONFIG} --sharedir OUTPUT_VARIABLE PG_CONFIG_SHAREDIR OUTPUT_STRIP_TRAILING_WHITESPACE)

find_package(BZip2 REQUIRED)
include(ExternalProject)

set(EXTERNAL_C_FLAGS "-g -Wall -O2 -fPIC -march=native -I${BIOSEQDB_DEPENDENCIES_PATH}/simde -DSIMDE_ENABLE_NATIVE_ALIASES")

ExternalProject_Add(bwamem
  SOURCE_DIR "${BIOSEQDB_DEPENDENCIES_PATH}/bwa"
  BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/bwamem"
  BUILD_COMMAND make CFLAGS=${EXTERNAL_C_FLAGS} CXXFLAGS=${EXTERNAL_C_FLAGS} libbwa.a
  BUILD_IN_SOURCE 0
  CONFIGURE_COMMAND cp -R ${BIOSEQDB_DEPENDENCIES_PATH}/bwa/. .
  UPDATE_COMMAND "" # Skip annoying updates for every build
  INSTALL_COMMAND ""
)

ExternalProject_Get_Property(bwamem install_dir)

add_library(libbwa SHARED IMPORTED)
set_target_properties(libbwa PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/bwamem/libbwa.a)
add_dependencies(libbwa bwamem)

file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/extension")

add_library(${POSTGRES_EXTENSION_NAME} SHARED
        ${SRC_FILES}
        ${CMAKE_CURRENT_SOURCE_DIR}/gen/config.h
        gen/extension.cpp
        )

add_dependencies(${POSTGRES_EXTENSION_NAME} libbwa)

target_include_directories(${POSTGRES_EXTENSION_NAME} PRIVATE ${PostgreSQL_TYPE_INCLUDE_DIR})
target_link_libraries(${POSTGRES_EXTENSION_NAME} PRIVATE ${PostgreSQL_LIBRARIES} ${BZIP2_LIBRARIES} libbwa)
target_link_libraries(${POSTGRES_EXTENSION_NAME} PRIVATE rpc)

install(TARGETS ${POSTGRES_EXTENSION_NAME} DESTINATION ${PG_CONFIG_PKGLIBDIR})
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/gen/${POSTGRES_EXTENSION_NAME}.control DESTINATION ${PG_CONFIG_SHAREDIR}/extension)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/gen/${POSTGRES_EXTENSION_NAME}--${POSTGRES_EXTENSION_VERSION}.sql DESTINATION ${PG_CONFIG_SHAREDIR}/extension)
