#include "startFunctions.cpp"

int main(){

    // Default functions
    // 2*x*exp(-x^2), analytical solution 1.0
    // 2*sin(1/x)^2, analytical solution M_PI

    std::cout << std::endl;

    // lower range, higher range, number of samples and analytical solution
    double parameters[] = {0.0, 1.0, 1000000.0, 0.0};

    // Define functions to be calculated
    // Note: the functions have to chosen such the f(x) > 0 for x > 0
    // and the sample number have to be even (Prerequisite for Simpson's method)
    std::string modified_gaussian_function = "4*sqrt(1-x^2)";
    parameters[3] = M_PI;
    obtainIntegralValue(&modified_gaussian_function, parameters);


    std::string oscillatory_sinus_function = "2*sin(1/x)^2";
    parameters[0] = 0.001;
    parameters[1] = 1000;
    parameters[3] = 3.14129;
    obtainIntegralValue(&oscillatory_sinus_function, parameters);


    return 0;


}
