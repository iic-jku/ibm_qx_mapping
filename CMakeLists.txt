cmake_minimum_required (VERSION 3.0)
project( ibm_qx_mapping )

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og -Wall -Wextra -Wpedantic -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -DNDEBUG -g0 -O3")

FIND_PACKAGE( Boost 1.50 COMPONENTS program_options REQUIRED )
file(GLOB_RECURSE SOURCES src/main.cpp src/util.cpp src/architecture_handling.cpp src/cost.cpp src/layer_handling.cpp src/circuit_property_handling.cpp src/node_handling.cpp src/mapping.cpp src/qasm/QASMparser.cpp src/qasm/QASMscanner.cpp src/qasm/QASMtoken.cpp)

add_executable(ibm_qx_mapping ${SOURCES})

include_directories(src src/qasm  ${Boost_INCLUDE_DIR})
target_link_libraries(ibm_qx_mapping boost_program_options)
