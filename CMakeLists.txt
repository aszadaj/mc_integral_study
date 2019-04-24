cmake_minimum_required(VERSION 3.13)
project(mc_integral_study)

set(CMAKE_CXX_STANDARD 14)

# This line is to find the NumPy version installed through pip
set(PYTHON_INCLUDE_DIR ${PYTHON_INCLUDE_DIR};/Library/Python/2.7/site-packages/numpy/core/include)

add_executable(mc_integral_study main.cpp fparser4/fparser.cc)

target_include_directories(mc_integral_study PRIVATE ${PYTHON_INCLUDE_DIR})
target_link_libraries(mc_integral_study ${PYTHON_LIBRARY})


# Debug for localizing which cmake variable to change.
macro(print_all_variables)
    message(STATUS "print_all_variables------------------------------------------{")
    get_cmake_property(_variableNames VARIABLES)
    foreach (_variableName ${_variableNames})
        message(STATUS "${_variableName}=${${_variableName}}")
    endforeach()
    message(STATUS "print_all_variables------------------------------------------}")
endmacro()

#print_all_variables()
