# Monte Carlo Integral Study

Project in the course SI2530 Computational Physics (KTH Royal Institute of Technology, 2019)


The software analyzes two integrals using Simpson's method, Simple Monte Carlo and a Metropolis algorithm. 
Main aspect is to investigate CPU time, that is how quickly the methods reaches an error level, the correlation time,
which sample is correlated with the distribution using the Metropolis method and comparisons between RMS and standard
error. The output is either a log on the console window or in form of a plot.

# How to run

The code is adapted to run with C++11 with ```matplotlibcpp``` library which have ```numpy``` as a dependency. The code
can be compiled with the input in terminal on macOS by

    gcc++ main.cpp
 
 when the compilation is ready, the executable is located in ```cmake-build-debug/mc_integral_study```. To run within
 the terminal, one can use the command
 
    /.cmake-build-debug/mc_integral_study
    
 Another way is to open the project, that is the main folder, that is ```mc_integral_study``` in CLion and to build
 and run with ```CMD+R```.
 

# Functions

The main class containing the relevant methods are in ```NumericalMethods.cpp```. Inside ```startFunctions.cpp``` are
four functions which uses ```NumericalMethods```-object, called ```calculations``` to analyze the integrals.

Obtains the integral value using N=10e6 samples and with random walk delta = 1.

    calculateIntegral(&calculations);

Analyzes errors by exporting a plot with RMS error and average standard error for various N.
   
    analyzeErrors(&calculations);

Obtains information in the console window how quickly the methods reaches an error level

    analyzeCPUTimes(&calculations);

Gets correlation times using only the Metropolis method and prints it to the console window. 

    analyzeCorrelationTime(&calculations);
    
