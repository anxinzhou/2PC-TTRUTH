cmake_minimum_required(VERSION 3.17)
project(ttruth)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE DEBUG)
add_subdirectory(ABY)
include_directories("./")

add_executable(ttruth main.cpp common/mpc_util.cpp mpc_ttruth/mpc_ttruth.cpp ttruth/ttruth.cpp ttruth/ttruth.h test.cpp test.h)

target_link_libraries(ttruth ABY::aby ENCRYPTO_utils::encrypto_utils)
