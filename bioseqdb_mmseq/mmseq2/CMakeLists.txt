cmake_minimum_required(VERSION 3.17)
project(mmseq2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(rpclib)
include_directories(rpclib/include)

find_package(Threads REQUIRED) # multithreading
find_package(PostgreSQL REQUIRED) # libpq - db connection

add_executable(mmseq2 main.cpp ../common/mmseq2lib.h dbconn.h dbconn.cpp mmseq2.h mmseq2.cpp)

target_link_libraries(mmseq2 PRIVATE PostgreSQL::PostgreSQL) 
target_link_libraries(mmseq2 PRIVATE Threads::Threads)
target_link_libraries(mmseq2 PRIVATE rpc)
