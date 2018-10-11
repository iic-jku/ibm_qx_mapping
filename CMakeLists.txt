cmake_minimum_required (VERSION 3.0)
project( ibm_qx_mapping )

SET(CMAKE_CXX_FLAGS "-Wall -g3 -O3 -D_REENTRANT -DVERSION=\"$(VERSION)\" -std=c++11")

add_executable(ibm_qx_mapping src/main.cpp)
include_directories(src)
