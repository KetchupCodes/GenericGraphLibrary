cmake_minimum_required(VERSION 3.4...3.18)
project(graphlibrary)

add_library(graphlib SHARED src/graphlib.cpp)
set_target_properties(graphlib PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_STANDARD 17)

add_subdirectory(pybind11)
pybind11_add_module(graphlibrary src/binder.cpp)
target_link_libraries(graphlibrary PUBLIC graphlib)
