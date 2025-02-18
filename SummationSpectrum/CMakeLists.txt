cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
PROJECT(main CXX C)
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src/)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
add_library(mylib SHARED main.cc ${sources} ${headers})
target_link_libraries(mylib ${ROOT_LIBRARIES})
add_executable(main main.cc ${sources} ${headers})
target_link_libraries(main mylib)