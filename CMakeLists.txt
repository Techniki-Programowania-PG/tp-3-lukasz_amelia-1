cmake_minimum_required(VERSION 3.14)
project(Sig)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(Sig Sig.cpp)

add_subdirectory(matplotplusplus)


target_link_libraries(Sig PUBLIC matplot)
