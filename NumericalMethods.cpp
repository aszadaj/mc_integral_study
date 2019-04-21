//
// Created by Antek Szadaj on 2019-02-19.
//

#include "NumericalMethods.h"

// Constructor
NumericalMethods::NumericalMethods(double parameters [], std::string * function) {

    // lower and higher limits
    a = parameters[0];
    b = parameters[1];

    samples = (int) parameters[2];
    standardDeviation = 0.0;
    analyticalSolution = parameters[3];
    timePassed = 0.0;

    // step size, h or abbreviated in the coursebook as \delta x
    h = (b-a)/samples;

    // Define which function is used from the string
    fparser.Parse(*function, "x");

}

// Integral calculation using the rectangular method (integral slicing)
void NumericalMethods::rectangular() {

    this->startClock();

    // calculate the integral value from i = 0 to i = n-1
    double x = a;

    while (x < b){

        integralValue += fparser.Eval(&x) * h;
        x += h;

    }

    this->stopClockAndPrintResults("Rectangular");

}

// Integral calculation using the trapezoidal method (integral slicing)
void NumericalMethods::trapezoidal() {

    this->startClock();

    // calculate the integral from i = 1 to i = n-1
    double x = a + h;

    while (x < b) {

        integralValue += fparser.Eval(&x) * h;
        x += h;

    }

    // add the boundary values of the trapezoidal method
    integralValue += 0.5 * (fparser.Eval(&a) + fparser.Eval(&b)) * h;

    this->stopClockAndPrintResults("Trapezoidal");

}

// Integral method using Simpson's method (integral slicing)
void NumericalMethods::simpson() {

    this->startClock();

    // calculate the integral from i = 1 to i = n-1
    double x = a + h;

    while (x < b) {

        // first line is either 4 or 2
        integralValue += 2*((int)((x-a)/h)%2+1) * fparser.Eval(&x) * 1/3.0 * h;
        x += h;

    }

    // add the boundary values of the Simpson's method
    integralValue += (fparser.Eval(&a) + fparser.Eval(&b)) * 1/3.0 * h;

    this->stopClockAndPrintResults("Simpson's");

}

// Integral calculation using Monte Carlo Mean Sample Value method
void NumericalMethods::MCMeanIntegral(){

    this->startClock();

    // randomized x-value and y-value for the function
    double x_rand, y_value;

    // standard deviation based on the randomized x-values
    double mean_value = 0.0;
    double squared_mean_value = 0.0;

    // Start seed and create distribution for the range [a,b]
    std::random_device rdX;
    std::mt19937 engX(rdX());
    std::uniform_real_distribution<> distrX(a,b);

    // Calculate the integral by evaluating f(x) at randomized x-value
    // and based on the randomized points, calculate the standard deviation

    for (int i = 0; i < samples; i++){

        x_rand = distrX(engX);
        y_value = fparser.Eval(&x_rand);
        integralValue += y_value * h;

        mean_value += y_value / samples;
        squared_mean_value += std::pow(y_value,2) / samples;

    }

    standardDeviation = std::sqrt(squared_mean_value - std::pow(mean_value,2));

    this->stopClockAndPrintResults("Monte Carlo mean sample");

}

// Integral calculation using Monte Carlo Hit or Miss method
void NumericalMethods::MCHitOrMissIntegral() {

    // Pairs of randomized points
    double x_rand, y_rand;

    // y-value for the function at x_rand
    double y_value;

    // standard deviation based on the randomized x-values
    double mean_value = 0.0;
    double squared_mean_value = 0.0;

    // Number of random samples which are below the function
    int samples_in_integral = 0;

    // Define the height of the rectangle from x-axis
    double y_max = this->getMaximumValue();

    // Start seed and distribution for x = [a,b] and y = [0.0, y_max], that is within the rectangle
    std::random_device rdX, rdY;
    std::mt19937 engX(rdX()), engY(rdY());
    std::uniform_real_distribution<> distrX(a, b), distrY(0.0, y_max);

    // Start at this point, to not include the time it takes to find the maximum point
    this->startClock();

    // Calculate how many randomized samples are within the area
    for (int i = 0; i < samples; i++){

        // Randomize a point
        x_rand = distrX(engX);
        y_rand = distrY(engY);

        y_value = fparser.Eval(&x_rand);

        // Based on the randomized points, calculate the variance
        mean_value += y_value / samples;
        squared_mean_value += std::pow(y_value, 2) / samples;

        // check if y-value is within the integral area
        if (y_rand < y_value)
            samples_in_integral++;
    }

    standardDeviation = std::sqrt(squared_mean_value - std::pow(mean_value,2));

    // the integral is the fraction of the hits under the function line
    integralValue = (b - a) * y_max * samples_in_integral/samples;

    this->stopClockAndPrintResults("Monte Carlo hit or miss");

}

// Calculate the height of the rectangle, by taking the maximum value of the integrand
double NumericalMethods::getMaximumValue(){

    double x = a;
    double y_max, y_value;

    while(x < b){

        y_value = fparser.Eval(&x);

        if (y_value > y_max)
            y_max = y_value;

        x += h;

    }

    return y_max;
}

// Calculate standard deviation for method which slices the integral
void NumericalMethods::calculateStandardDeviation(){

    if (standardDeviation != 0.0)
        return;

    double y_value, mean_value, mean_value_squared;
    double x_value = a;

    for (int i = 0; i < samples; i++){

        y_value = fparser.Eval(&x_value);
        mean_value += y_value / samples;
        mean_value_squared += std::pow(y_value,2) / samples;
        x_value += h;

    }

    standardDeviation = std::sqrt(mean_value_squared - std::pow(mean_value,2));

}

// Start the time and clear the integral value
void NumericalMethods::startClock() {

    timePassed = (double) clock();
    integralValue = 0.0;

}

// Stop clock and print the results
void NumericalMethods::stopClockAndPrintResults(std::string methodName) {

    timePassed = ((double) clock() - timePassed) * 1000.0 / CLOCKS_PER_SEC;

    this->calculateStandardDeviation();

    std::cout << methodName << " method:" << std::endl;
    std::cout << "Integral value: \t\t" << integralValue << std::endl;
    std::cout << "Error: \t\t\t\t\t" << std::abs(integralValue - analyticalSolution) << ", ";
    std::cout << std::abs(integralValue - analyticalSolution)*100/analyticalSolution << " %" << std::endl;
    std::cout << "Time: \t\t\t\t\t" << timePassed << " ms" << std::endl;

    // Note: The standard error is defined as \sigma / \sqrt{N}, with N being number of samples
    std::cout << "Standard error: \t\t" << standardDeviation/std::sqrt(samples) << ", ";
    std::cout << standardDeviation/std::sqrt(samples)*100/analyticalSolution << " %" << std::endl << std::endl;
}

double NumericalMethods::getA() const {
    return a;
}

double NumericalMethods::getB() const {
    return b;
}

int NumericalMethods::getSamples() const {
    return samples;
}

double NumericalMethods::getAnalyticalSolution() const {
    return analyticalSolution;
}
