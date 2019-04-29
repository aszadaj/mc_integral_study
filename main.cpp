#include "startFunctions.cpp"

int main() {

    double long analyticalSolution;

    float lowerRange, higherRange;
    float *ranges[] = {&lowerRange, &higherRange};

    std::string integrandFunctionString, PDFFunctionString;
    std::string *functionStrings[] = {&integrandFunctionString, &PDFFunctionString};

    // the functions are such that the PDF have to be normalized in order to calculate the integrals
    integrandFunctionString = "x*exp(-x^2)";
    PDFFunctionString = "2*exp(-x^2)/sqrt(3.1415926536)"; //"2*exp(-x^2)/sqrt(3.1415926536)"
    lowerRange = 0.0;
    higherRange = 10000.0;
    analyticalSolution = 0.5;

    obtainIntegralValue(&analyticalSolution, ranges, functionStrings);

    return 0;

}