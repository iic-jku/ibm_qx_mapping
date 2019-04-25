cmake_minimum_required (VERSION 3.0)
project( ibm_qx_mapping )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

SET(CMAKE_CXX_FLAGS "-Wall -g -O2 -DNDEBUG -DVERSION=\"$(VERSION)\" -std=c++11")


add_executable(ibm_qx_mapping src/main.cpp src/QASMparser.cpp src/QASMscanner.cpp src/QASMtoken.cpp)
include_directories(src)
