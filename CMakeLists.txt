cmake_minimum_required(VERSION 3.22)
project(Homework3_int2048)

set(CMAKE_CXX_STANDARD 14)
set(CMake_CXX_FLAGS "-g")

include_directories(.)

add_executable(Homework3_int2048
        int2048.cpp
        int2048.h)
