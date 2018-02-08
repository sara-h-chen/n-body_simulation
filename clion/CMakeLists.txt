cmake_minimum_required(VERSION 3.6)
project(clion)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp)
add_executable(clion ${SOURCE_FILES})

#find_package(PythonLibs 2.7)
#target_include_directories(clion PRIVATE ${PYTHON_INCLUDE_DIRS})
#target_link_libraries(clion ${PYTHON_LIBRARIES})