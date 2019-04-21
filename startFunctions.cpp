//
// Created by Antek Szadaj on 2019-03-13.
//

#include "startFunctions.h"
#include <string>


void obtainIntegralValue(std::string * function, double parameters []){

    NumericalMethods calculations = NumericalMethods(parameters, &*function);

    std::cout << "Integral for x = [" << calculations.getA() << "," << calculations.getB() << "] ";
    std::cout << "of f(x) = " << *function << std::endl;
    std::cout << "Analytical solution: " << calculations.getAnalyticalSolution() << std::endl;
    std::cout << "Number of samples: " << calculations.getSamples() << std::endl << std::endl;


    // Integral slicing methods
    calculations.rectangular();
    calculations.trapezoidal();
    calculations.simpson();


    // Probability integrals
    calculations.MCMeanIntegral();
    calculations.MCHitOrMissIntegral();

    std::cout << std::endl << std::endl;


}