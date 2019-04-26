//
// Created by Antek Szadaj on 2019-02-19.
//

#ifndef REPORT_CODE_NUMERICALMETHODS_H
#define REPORT_CODE_NUMERICALMETHODS_H

#include <iostream>
#include "fparser4/fparser.hh"
#include <random>


class NumericalMethods {

    double stepSize;

    FunctionParser integrandFunction;
    FunctionParser densityFunction;
    unsigned int samples;

public:



    double lowerLimit;
    double higherLimit;
    unsigned int sampleLevel;
    double analyticalSolution;
    double integralResult;



    double standardError;
    double errorFromRealResult;
    double correlationTime;
    double timePassed;
    double errorLevel;
    float delta;

    double randomWalkTime;

    bool printMessage;
    bool CPUTimeAnalysis;

    NumericalMethods(double * [], std::string * []);

    void simpson();
    void simpleMonteCarlo();
    void metropolis();

    void startClock(double = 0.0);
    void stopClock();
    void printResults(std::string);

    double getArgMax();
    void getStandardError(double [], bool = true);
    double getAutocorrelationValue(double [], int);
    void createRandomWalk(double []);
    unsigned int determineCorrelatedStep(double []);
    unsigned int getSamples() const;
    void setSamples(unsigned int samples);

    void resetValues();

};


#endif //REPORT_CODE_NUMERICALMETHODS_H
