#include "startFunctions.cpp"

int main() {

    double lowerRange, higherRange, numberOfSamples, analyticalSolution;
    double *parameters[] = {&lowerRange, &higherRange, &numberOfSamples, &analyticalSolution};

    std::string mainFunctionString, weightFunctionString, densityFunctionString;
    std::string *functionStrings[] = {&mainFunctionString, &weightFunctionString, &densityFunctionString};

    mainFunctionString = "x";
    weightFunctionString = "exp(-x)";
    densityFunctionString = "0.5*exp(-x/2)";

    lowerRange = 0.0;
    higherRange = 10000.0;
    analyticalSolution = 1.0;
    numberOfSamples = 1.0e5;

    obtainIntegralValue(parameters, functionStrings);

    return 0;

}