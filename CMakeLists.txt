cmake_minimum_required (VERSION 3.0)
project( ibm_qx_mapping )

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og -Wall -Wextra -Wpedantic -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -DNDEBUG -g0 -O3")

file(GLOB_RECURSE SOURCES src/main.cpp src/QASMparser.cpp src/QASMscanner.cpp src/QASMtoken.cpp)
add_executable(ibm_qx_mapping ${SOURCES})
include_directories(src)
