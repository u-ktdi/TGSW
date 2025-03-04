cmake_minimum_required(VERSION 3.10)
project(TGSWProject)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(tgsw SHARED tgsw/tgsw.cpp)
add_executable(test_tgsw tgsw/test_tgsw.cpp)
target_link_libraries(test_tgsw PRIVATE tgsw)