# Monte Carlo Integral Study

Project in the course SI2530 Computational Physics (KTH Royal Institute of Technology)


The code is built to analyze two integrals using Simpson's method, Simple Monte Carlo and a Metropolis algorithm.
Investigation such as error level, correlation time (for Metropolis method), RMS and standard errors and general 
comparisons of those methods. 

# How to run

The code is adapted to run with C++11 with ```matplotlibcpp``` library which have ```numpy``` as a dependency. The code
can be compiled with the input in terminal on macOS by

    gcc++ main.cpp
 
 when the compilation is ready, the executable is located in ```cmake-build-debug/mc_integral_study```. To run within
 the terminal, one can use the command
 
    /.cmake-build-debug/mc_integral_study
    
 Another way is to open the project, that is the main folder, that is ```mc_integral_study``` in CLion and to build
 and run with ```CMD+R```.


# Output

The output is in form of console log or a plot using ```matplotlibcpp```, which are used as a reference in the report
(to be appeared).