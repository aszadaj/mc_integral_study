#include "startFunctions.cpp"

int main() {

    double  analyticalSolution, delta;
    double  *parameters[] = {&analyticalSolution, &delta};

    std::string integrandFunctionString, densityFunctionString;
    std::string *functionStrings[] = {&integrandFunctionString, &densityFunctionString};

    integrandFunctionString = "exp(-x^2)"; // x*exp(-x^2)
    densityFunctionString = "exp(-x)*exp(1)/(exp(1)-1)"; // 2*exp(-x^2)/sqrt(3.14159)

    // analytical solution for
    analyticalSolution = 0.74689;
    delta = 1;

    obtainIntegralValue(parameters, functionStrings);

    return 0;

}