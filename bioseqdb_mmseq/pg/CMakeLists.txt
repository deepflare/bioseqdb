
cmake_minimum_required(VERSION 3.17)
project(mmseq2 LANGUAGES CXX)

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

include_directories(${CMAKE_CURRENT_BINARY_DIR})
configure_file(config.h.in config.h @ONLY)
add_library(mmseq2 SHARED extension.cpp config.h ../common/mmseq2lib.h)

# Add RPC support
add_subdirectory(rpclib)
include_directories(rpclib/include)
target_link_libraries(mmseq2 PRIVATE rpc)

# Add Postgres support
find_package(PostgreSQL REQUIRED)
target_link_libraries(mmseq2 PRIVATE PostgreSQL::PostgreSQL)

find_program(PG_CONFIG pg_config)
execute_process(COMMAND ${PG_CONFIG} --pkglibdir OUTPUT_VARIABLE PG_CONFIG_PKGLIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${PG_CONFIG} --sharedir OUTPUT_VARIABLE PG_CONFIG_SHAREDIR OUTPUT_STRIP_TRAILING_WHITESPACE)

# Actually install files
install(TARGETS mmseq2 DESTINATION ${PG_CONFIG_PKGLIBDIR})
install(FILES mmseq2.control DESTINATION ${PG_CONFIG_SHAREDIR}/extension)
install(FILES mmseq2--0.0.0.sql DESTINATION ${PG_CONFIG_SHAREDIR}/extension)
