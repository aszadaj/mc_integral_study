//
// Created by Antek Szadaj on 2019-03-13.
//

#include "startFunctions.h"
#include <string>


void obtainIntegralValue(double * parameters [], std::string * functionStrings []){

    NumericalMethods calculations = NumericalMethods(parameters, functionStrings);

    std::cout << std::endl;

    std::cout << "Integral: int(" << functionStrings[0] << ","<< calculations.getA() <<","<< calculations.getB();
    std::cout << ",x) / int("<< *functionStrings[1] <<","<< calculations.getA() <<","<< calculations.getB() << ",x)";
    std::cout << "Density function: " << *functionStrings[2] << std::endl;
    std::cout << "Analytical solution: " << calculations.getAnalyticalSolution() << std::endl;
    std::cout << "Number of samples: " << calculations.getSamples() << std::endl << std::endl;

    // Integral slicing methods
    calculations.rectangularMethod();
    calculations.trapezoidalMethod();
    calculations.simpsonMethod();

    // Probability integrals
    calculations.monteCarloMeanSampleMethod();
    calculations.monteCarloHitOrMissMethod();
    calculations.metropolisMethod();

    std::cout << std::endl << std::endl;


}