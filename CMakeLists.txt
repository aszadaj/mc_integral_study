cmake_minimum_required(VERSION 3.13)
project(report_code)

set(CMAKE_CXX_STANDARD 14)

add_executable(mc_integral_study main.cpp fparser4/fparser.cc)