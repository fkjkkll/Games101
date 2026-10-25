cmake_minimum_required(VERSION 3.10)
project(BezierCurve)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenCV REQUIRED)

add_executable(BezierCurve main.cpp)
if(MSVC)
    target_compile_options(BezierCurve PRIVATE /utf-8)
endif()
target_link_libraries(BezierCurve ${OpenCV_LIBRARIES})
