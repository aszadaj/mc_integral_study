//
// Created by Antek Szadaj on 2019-02-19.
//

#ifndef REPORT_CODE_NUMERICALMETHODS_H
#define REPORT_CODE_NUMERICALMETHODS_H

#include <iostream>
#include "fparser4/fparser.hh"
#include <random>
#include "matplotlibcpp.h"


class NumericalMethods {

public:

    double a;
    double b;
    unsigned int samples;
    double analyticalSolution;

private:

    double stepSize;
    double timePassed;

    FunctionParser mainFunction;
    FunctionParser weightFunction;
    FunctionParser densityFunction;

    double integralResult;

public:


    NumericalMethods(double * [], std::string * []);

    void simpson();
    void simpleMonteCarlo();
    void metropolis();

    void startClock();
    void stopClockAndPrintResults(std::string, double = 0.0);

    double getXPositionForMaximumValue();
    double calculateStandardError(double []);


};


#endif //REPORT_CODE_NUMERICALMETHODS_H
