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
    double getA() const;

    double getB() const;

    unsigned int getSamples() const;

    double getAnalyticalSolution() const;

private:

    // h = separation on x-axis
    double h;
    double timePassed;

    // Ranges for the integral and its sum
    // default ranges 0.00001, 10000

    FunctionParser numeratorFunction;
    FunctionParser denominatorFunction;
    FunctionParser densityFunction;

    double a;
    double b;

    unsigned int samples;
    double standardDeviationNumerator;
    double standardDeviationDenominator;
    double analyticalSolution;

    double integralResult;



public:


    NumericalMethods(double * [], std::string * []);

    void rectangularMethod();
    void trapezoidalMethod();
    void simpsonMethod();

    void monteCarloMeanSampleMethod();
    void monteCarloHitOrMissMethod();
    void metropolisMethod();

    void startClock();
    void stopClockAndPrintResults(std::string);

    double getMaximumValueIntegrand(bool = true);
    double getXPositionForMaximumValue();
    void calculateStandardDeviation();


};


#endif //REPORT_CODE_NUMERICALMETHODS_H
