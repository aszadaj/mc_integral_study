#include "startFunctions.cpp"

int main() {

    double long analyticalSolution;

    // Analyze int(x*exp(-x^2),0,inf) with density function
    // p(x) = 2*exp(-x^2) / sqrt(pi)
    analyticalSolution = 0.5;
//    obtainIntegralValue(&analyticalSolution);

    // Analyze int(sin(1/x),0,inf) with density function
    // 1. p_n(x) = 0.5 * (n * pi^2 (n + 1)) * sin ( n * pi^2 (n + 1) x - n * pi )
    // where x is on the range [1 / ( (n+1)pi ), 1 / (n * pi)] with integer n > 1
    // 2. p(x) = pi / (1 + pi) * x * exp(1/pi - x)
    // for x > 1/pi
    analyticalSolution = M_PI/2;
    obtainIntegralValue(&analyticalSolution, false);

    return 0;

}