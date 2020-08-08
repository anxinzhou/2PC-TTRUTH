cmake_minimum_required(VERSION 3.17)
project(ttruth)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(ABY)
include_directories("./")

add_executable(ttruth main.cpp common/mpc_util.cpp mpc_ttruth/mpc_ttruth.cpp)

target_link_libraries(ttruth ABY::aby ENCRYPTO_utils::encrypto_utils)
