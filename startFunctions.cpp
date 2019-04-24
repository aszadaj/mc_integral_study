//
// Created by Antek Szadaj on 2019-03-13.
//

#include "startFunctions.h"
#include <string>


void obtainIntegralValue(double * parameters [], std::string * functionStrings []){

    NumericalMethods calculations = NumericalMethods(parameters, functionStrings);

    printIntegralInformation(functionStrings, &calculations);

    calculations.simpson();
    calculations.simpleMonteCarlo();
    calculations.metropolis();


}

void printIntegralInformation(std::string * functionStrings [], NumericalMethods * calculations){

    std::cout << std::endl << "Integral: int(" << *functionStrings[0] << "*"<< *functionStrings[1] <<",";
    std::cout << calculations->a << "," << calculations->b << ",x) / int("<< *functionStrings[1] <<",";
    std::cout << calculations->a <<","<< calculations->b << ",x)" << std::endl << "Density function: ";
    std::cout << *functionStrings[2] << std::endl << "Analytical solution: " << calculations->analyticalSolution;
    std::cout << std::endl << "Number of samples: " << calculations->samples << std::endl << std::endl;

}