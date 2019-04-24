#include "startFunctions.cpp"

int main(){

    double lowerRange, higherRange, numberOfSamples, analyticalSolution;
    double * parameters [] = {&lowerRange, &higherRange, &numberOfSamples, &analyticalSolution};

    std::string numeratorFunctionString, denominatorFunctionString, densityFunctionString;
    std::string * functionStrings [] = {&numeratorFunctionString, &denominatorFunctionString, &densityFunctionString};

    numeratorFunctionString = "x*exp(-x)";
    denominatorFunctionString = "exp(-x)";
    densityFunctionString = "0.5*exp(-x/2)";
    lowerRange = 0.0;
    higherRange = 10000.0;
    analyticalSolution = 1.0;
    numberOfSamples = 1.0e6;

    obtainIntegralValue(parameters, functionStrings);

    return 0;

}



// P(x) = int(p(y),-inf, x, y)

//-exp(-x/2) {0, x} = 1-exp(-x/2) = r -> ln (1-r) = -x/2 -> x = -2*ln(r)