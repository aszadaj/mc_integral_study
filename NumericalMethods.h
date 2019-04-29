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

    double lowerLimit, higherLimit, stepSize;

    unsigned int samples;

    FunctionParser integrandFunction;
    FunctionParser PDFFunction;

public:

    double long analyticalSolution;
    double integralResult, standardError, errorFromRealResult, correlationTime, integralTime, errorLevel, randomWalkTime;
    unsigned int sampleLevel, rejectedSamples;
    float delta;

    // maybe delete randomWalkTime


    bool printMessage, CPUTimeAnalysis;


    NumericalMethods(const double long *, float * [], std::string * []);

    void simpson();
    void simpleMonteCarlo();
    void metropolis();

    void startClock(double = 0.0);
    void stopClock();
    void printResults(std::string);

    unsigned int determineCorrelatedStep(double []);

    void resetValues();

    void setSamples(unsigned int samples);
    void setLowerLimit(double lowerLimit);
    void setHigherLimit(double higherLimit);

    double getLowerLimit() const;
    double getHigherLimit() const;
    double getArgMax();
    void getStandardError(double []);
    double getAutocorrelationValue(double [], int);
    unsigned int getSamples() const;


};


#endif //REPORT_CODE_NUMERICALMETHODS_H
