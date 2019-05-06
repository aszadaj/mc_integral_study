#include "startFunctions.cpp"

int main() {

    double long analyticalSolution;

    // Analyze int(x*exp(-x^2),0,inf) with density function
    // p(x) = 2*exp(-x^2) / sqrt(pi)

    analyticalSolution = 0.5;
    obtainIntegralValue(&analyticalSolution);

    // Analyze int(sin^2(1/x),0,inf) with density function
    // 1. p_n(x) = N * sin ( n * pi^2 (n + 1) x - n * pi )
    // where x is on the range [1 / ( (n+1)pi ), 1 / (n * pi)] with integer n > 0
    // 2. p(x) = N * (pi * x - 1) for x on [1/pi, 2/pi]
    // 3. p(x) = N * exp(2/pi - x) for x > 2/pi

    analyticalSolution = M_PI/2;
    obtainIntegralValue(&analyticalSolution, false);

    return 0;

}