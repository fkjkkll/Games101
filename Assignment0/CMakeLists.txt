cmake_minimum_required(VERSION 3.10)
project(Transformation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)
get_target_property(_EIGEN3_INC Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
include_directories(${_EIGEN3_INC}/..)

add_executable(Transformation main.cpp)
if(MSVC)
    target_compile_options(Transformation PRIVATE /utf-8)
endif()
