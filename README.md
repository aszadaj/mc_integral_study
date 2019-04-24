#Monte Carlo Integral Study

Project in the course SI2530 Computational Physics (KTH Royal Institute of Technology)


The code is built to analyze the methods described in the ```project_formulation.pdf```.

# How to run

The code is adapted to run with C++11 and the compilation can either be done with GNU compiler or other. On macOS console (terminal) the compilation is done bu

    gcc++ main.cpp
 
 when the compilation is ready, the executable is located in ```cmake-build-debug/mc_integral_study```. To run within
 the terminal, one can use the command
 
    /.cmake-build-debug/mc_integral_study
    
 Another way is to open the project, that is the main folder, that is ```mc_integral_study``` in CLion and to build
 and run with ```CMD+R```.


# Output

The output gives various information about the integrated methods, which are Rectangular, Trapezoidal, Simpson's, Monte
Carlo mean sample method and Monte Carlo hit or miss method.