cmake_minimum_required (VERSION 3.0)
project( ibm_qx_mapping )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

SET(CMAKE_CXX_FLAGS "-Wall -g3 -O3 -D_REENTRANT -DVERSION=\"$(VERSION)\" -std=c++11")

FIND_PACKAGE(MPFR REQUIRED)
FIND_PACKAGE(MPFR++ REQUIRED)
INCLUDE_DIRECTORIES(${MPFR_INCLUDE_DIRS})
SET(LIBS ${LIBS} ${MPFR_LIBRARIES})


add_executable(ibm_qx_mapping src/main.cpp src/QASMparser.cpp src/QASMscanner.cpp src/QASMtoken.cpp)
include_directories(src)
target_link_libraries(ibm_qx_mapping mpfr)
