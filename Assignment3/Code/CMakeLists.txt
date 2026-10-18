cmake_minimum_required(VERSION 3.10)
project(Rasterizer)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenCV REQUIRED)
find_package(Eigen3 REQUIRED)
get_target_property(_EIGEN3_INC Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
include_directories(${_EIGEN3_INC}/..)

include_directories(./include)

add_executable(Rasterizer main.cpp rasterizer.hpp rasterizer.cpp global.hpp Triangle.hpp Triangle.cpp Texture.hpp Texture.cpp Shader.hpp OBJ_Loader.h)
if(MSVC)
    target_compile_options(Rasterizer PRIVATE /utf-8)
endif()
target_compile_definitions(Rasterizer PRIVATE RASTERIZER_MODEL_DIR="${CMAKE_CURRENT_SOURCE_DIR}/models")
target_link_libraries(Rasterizer ${OpenCV_LIBRARIES})
