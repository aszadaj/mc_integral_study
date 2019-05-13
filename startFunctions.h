//
// Created by Antek Szadaj on 2019-03-13.
//

#ifndef REPORT_CODE_STARTFUNCTIONS_H
#define REPORT_CODE_STARTFUNCTIONS_H

#include "NumericalMethods.cpp"
#include "matplotlibcpp.h"

void obtainIntegralValue(double long *, bool = true);

void calculateIntegrals(NumericalMethods *);

void analyzeErrors(NumericalMethods *);
void analyzeCorrelationTime(NumericalMethods *);
void analyzeCPUTimes(NumericalMethods *);

void exportErrorPlot (std::vector<double> *, std::vector<std::vector <double>> *, int, float, bool, int);
void printIntegralInformation(NumericalMethods *, bool);

#endif //REPORT_CODE_STARTFUNCTIONS_H