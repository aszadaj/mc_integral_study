//
// Created by Antek Szadaj on 2019-02-19.
//

#ifndef REPORT_CODE_NUMERICALMETHODS_H
#define REPORT_CODE_NUMERICALMETHODS_H

#include <iostream>
#include "fparser4/fparser.hh"
#include <random>


class NumericalMethods {

public:
    double getA() const;

    double getB() const;

    int getSamples() const;

    double getAnalyticalSolution() const;

private:

    // h = separation on x-axis
    double h;
    double timePassed;

    // Ranges for the integral and its sum
    // default ranges 0.00001, 10000

    FunctionParser fparser;
    double a;
    double b;

    int samples;
    double standardDeviation;
    double analyticalSolution;

    double integralValue;



public:


    NumericalMethods(double [], std::string *);

    void rectangular();
    void trapezoidal();
    void simpson();

    void MCMeanIntegral();
    void MCHitOrMissIntegral();

    void startClock();
    void stopClockAndPrintResults(std::string);

    double getMaximumValue();
    void calculateStandardDeviation();


};


#endif //REPORT_CODE_NUMERICALMETHODS_H
