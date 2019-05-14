
#ifndef REPORT_CODE_NUMERICALMETHODS_H
#define REPORT_CODE_NUMERICALMETHODS_H

#include <iostream>
#include <random>
#include <sstream>
#include "matplotlibcpp.h"


class NumericalMethods {

    unsigned int samples;
    double lowerLimit, higherLimit, stepSize;

public:

    double long analyticalSolution;
    double integralResult, standardError, errorFromRealResult, integralTime, errorLevel;
    unsigned int rejectedSamples;
    float delta;

    bool printMessage, CPUTimeAnalysis, simpleIntegral, correlationTimeAnalysis;

    explicit NumericalMethods(bool = true);

    void simpson();
    void simpleMonteCarlo();
    void metropolis();

    void startClock(double = 0.0);
    void stopClock();
    void printResults(std::string);

    void exportCorrelationTimePlot(double []);
    void createRandomWalk(double []);

    void resetValues();

    void setSamples(unsigned int samples);

    double getMainFunction(double *);
    double getPDF(double *);
    double getSampledPDFValue(double *, double *, double *);

    double getLowerLimit() const;
    double getHigherLimit() const;
    double getRandomWalkStartValue();
    void getStandardError(double []);
    double getAutocorrelationValue(double [], int);
    unsigned int getSamples() const;

    std::string getFunctionName();

};


#endif //REPORT_CODE_NUMERICALMETHODS_H
