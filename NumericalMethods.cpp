//
// Created by Antek Szadaj on 2019-02-19.
//

#include "NumericalMethods.h"

namespace plt = matplotlibcpp;

// Constructor
NumericalMethods::NumericalMethods(double * parameters [], std::string * stringFunctions []) {

    // lower and higher limits
    a = *parameters[0];
    b = *parameters[1];

    samples = (int) *parameters[2];
    analyticalSolution = *parameters[3];

    standardDeviationNumerator = 0.0;
    standardDeviationDenominator = 0.0;
    timePassed = 0.0;
    integralResult = 0.0;

    // step size h, for slicing integrals
    h = (b-a)/samples;

    // Define numerator and denominator functions
    numeratorFunction.Parse(*stringFunctions[0], "x");
    denominatorFunction.Parse(*stringFunctions[1], "x");
    densityFunction.Parse(*stringFunctions[2], "x");

}


// Integral calculation using the rectangular method (integral slicing)
void NumericalMethods::rectangularMethod() {

    this->startClock();

    double numeratorResult = 0.0;
    double denominatorResult = 0.0;

    // calculate the integral value from i = 0 to i = n-1
    double x = a;

    while (x < b){

        numeratorResult += numeratorFunction.Eval(&x) * h;
        denominatorResult += denominatorFunction.Eval(&x) * h;
        x += h;

    }

    integralResult = numeratorResult/denominatorResult;

    this->stopClockAndPrintResults("Rectangular");

}


// Integral calculation using the trapezoidal method (integral slicing)
void NumericalMethods::trapezoidalMethod() {

    this->startClock();

    double numeratorResult = 0.0;
    double denominatorResult = 0.0;

    // calculate the integral from i = 1 to i = n-1
    double x = a + h;

    while (x < b) {

        numeratorResult += numeratorFunction.Eval(&x) * h;
        denominatorResult += denominatorFunction.Eval(&x) * h;
        x += h;

    }

    // add the boundary values of the trapezoidal method
    numeratorResult += 0.5 * (numeratorFunction.Eval(&a) + numeratorFunction.Eval(&b)) * h;
    denominatorResult += 0.5 * (denominatorFunction.Eval(&a) + denominatorFunction.Eval(&b)) * h;

    integralResult = numeratorResult/denominatorResult;

    this->stopClockAndPrintResults("Trapezoidal");

}


// Integral method using Simpson's method (integral slicing)
void NumericalMethods::simpsonMethod() {

    this->startClock();

    double numeratorResult = 0.0;
    double denominatorResult = 0.0;

    // calculate the integral from i = 1 to i = n-1
    double x = a + h;

    while (x < b) {

        // first line is either 4 or 2
        numeratorResult += 2*((int)((x-a)/h)%2+1) * numeratorFunction.Eval(&x) * 1/3.0 * h;
        denominatorResult += 2*((int)((x-a)/h)%2+1) * denominatorFunction.Eval(&x) * 1/3.0 * h;
        x += h;

    }

    // add the boundary values of the Simpson's method
    numeratorResult += (numeratorFunction.Eval(&a) + numeratorFunction.Eval(&b)) * 1/3.0 * h;
    denominatorResult += (denominatorFunction.Eval(&a) + denominatorFunction.Eval(&b)) * 1/3.0 * h;

    integralResult = numeratorResult/denominatorResult;

    this->stopClockAndPrintResults("Simpson's");

}


// Integral calculation using Monte Carlo Mean Sample Value method
void NumericalMethods::monteCarloMeanSampleMethod(){

    this->startClock();

    // randomized x-value
    double x_rand;

    // standard deviation based on the randomized x-values
    double mean_value_numerator = 0.0;
    double squared_mean_value_numerator = 0.0;

    double mean_value_denominator = 0.0;
    double squared_mean_value_denominator = 0.0;

    double standardDeviationDenominator;

    double integralNumerator = 0.0;
    double integralDenominator = 0.0;

    // Start seed and create distribution for the range [0,1]
    std::random_device randomDevice;
    std::mt19937 randomEngine(randomDevice());
    std::uniform_real_distribution<> distribution(0.0, 1.0);

    // Calculate the integral by evaluating f(x) at randomized x-value
    // and based on the randomized points, calculate the standard deviation
    // additionally, use importance sampling, to reduce the variance

    for (int i = 0; i < samples; i++){

        // x randomized according to the pdf = exp(-x) from 0 to infinity, x = -ln(r)
        x_rand = -2.0*log(distribution(randomEngine));

        integralNumerator += (numeratorFunction.Eval(&x_rand)/densityFunction.Eval(&x_rand))/samples;
        integralDenominator += (denominatorFunction.Eval(&x_rand)/densityFunction.Eval(&x_rand))/samples;

        mean_value_numerator += numeratorFunction.Eval(&x_rand)/samples;
        squared_mean_value_numerator += std::pow(numeratorFunction.Eval(&x_rand),2)/samples;

        mean_value_denominator += denominatorFunction.Eval(&x_rand)/samples;
        squared_mean_value_denominator += std::pow(denominatorFunction.Eval(&x_rand),2)/samples;

    }

    integralResult = (integralNumerator/integralDenominator);

    standardDeviationNumerator = std::sqrt(squared_mean_value_numerator - std::pow(mean_value_numerator,2));
    standardDeviationDenominator = std::sqrt(squared_mean_value_denominator - std::pow(mean_value_denominator,2));

    this->stopClockAndPrintResults("Monte Carlo mean sample");

}


// Integral calculation using Monte Carlo Hit or Miss method
void NumericalMethods::monteCarloHitOrMissMethod() {

    // Pairs of randomized points
    double x_rand_numerator, y_rand_numerator;
    double x_rand_denominator, y_rand_denominator;

    double integralResultNumerator, integralResultDenominator;

    // y-value for the function at x_rand
    double y_value_numerator, y_value_denominator;

    // standard deviation based on the randomized x-values
    double mean_value_numerator = 0.0;
    double squared_mean_value_numerator = 0.0;

    double mean_value_denominator = 0.0;
    double squared_mean_value_denominator = 0.0;

    // Number of random samples which are below the function
    unsigned int samples_in_integral_numerator = 0;
    unsigned int samples_in_integral_denominator = 0;

    // Define the height of the rectangle from x-axis
    double y_max_numerator = this->getMaximumValueIntegrand();
    double y_max_denominator = this->getMaximumValueIntegrand(false);

    // Start seed and distribution for x = [a,b] and y = [0.0, y_max], that is within the rectangle
    std::random_device randomDeviceXNumerator, randomDeviceYNumerator;
    std::mt19937 randomEngineXNumerator(randomDeviceXNumerator()), randomEngineYNumerator(randomDeviceYNumerator());
    std::uniform_real_distribution<> distributionXNumerator(a, b), distributionYNumerator(0.0, y_max_numerator);


    std::random_device randomDeviceXDenominator, randomDeviceYDenominator;
    std::mt19937 randomEngineXDenominator(randomDeviceXDenominator()), randomEngineYDenominator(randomDeviceYDenominator());
    std::uniform_real_distribution<> distributionXDenominator(a, b), distributionYDenominator(0.0, y_max_denominator);

    // Start at this point, to not include the time it takes to find the maximum point
    this->startClock();

    // Calculate how many randomized samples are within the area
    for (unsigned int i = 0; i < samples; i++){

        // Randomize a point
        x_rand_numerator = distributionXNumerator(randomEngineXNumerator);
        y_rand_numerator = distributionYNumerator(randomEngineYNumerator);

        x_rand_denominator = distributionXDenominator(randomEngineXDenominator);
        y_rand_denominator = distributionYDenominator(randomEngineYDenominator);

        y_value_numerator = numeratorFunction.Eval(&x_rand_numerator);
        y_value_denominator = denominatorFunction.Eval(&x_rand_denominator);

        // Based on the randomized points, calculate the variance
        mean_value_numerator += y_value_numerator / samples;
        squared_mean_value_numerator += std::pow(y_value_numerator, 2) / samples;

        // Based on the randomized points, calculate the variance
        mean_value_denominator += y_value_denominator / samples;
        squared_mean_value_denominator += std::pow(y_value_denominator, 2) / samples;

        // check if y-value is within the integral area
        if (y_rand_numerator < y_value_numerator)
            samples_in_integral_numerator++;

        if (y_rand_denominator < y_value_denominator)
            samples_in_integral_denominator++;

    }

    standardDeviationNumerator = std::sqrt(squared_mean_value_numerator - std::pow(mean_value_numerator,2));
    standardDeviationDenominator = std::sqrt(squared_mean_value_denominator - std::pow(mean_value_denominator,2));

    // the integral is the fraction of the hits under the function line
    integralResultNumerator = (b - a) * y_max_numerator * samples_in_integral_numerator/samples;
    integralResultDenominator = (b - a) * y_max_denominator * samples_in_integral_denominator/samples;

    integralResult = integralResultNumerator/integralResultDenominator;

    this->stopClockAndPrintResults("Monte Carlo hit or miss");

}

// Generate samples using Metropolis algorithm
void NumericalMethods::metropolisMethod() {

    double * x_i = new double[samples];

    double x_trial, delta, delta_i, transitionProb;

    double integralNumerator = 0.0;
    double integralDenominator = 0.0;

    double mean_value_numerator = 0.0;
    double squared_mean_value_numerator = 0.0;

    double mean_value_denominator = 0.0;
    double squared_mean_value_denominator = 0.0;

    double standardDeviationDenominator;

    // the size of the random walk
    delta = 1.5;

    // number of samples which have been accepted by the algorithm
    unsigned int acceptedSamples = 0;

    // Start seed and create distribution for the range [0,1]
    std::random_device randomDevice;
    std::mt19937 randomEngine(randomDevice());
    std::uniform_real_distribution <> distribution(0, 1.0);

    // Start the
    x_i[0] = this->getXPositionForMaximumValue();

    this-> startClock();

    for (unsigned int i = 0; i < samples; i++){

        delta_i = (2.0 * distribution(randomEngine) - 1.0) * delta;
        x_trial = x_i[i] + delta_i;
        transitionProb = densityFunction.Eval(&x_trial)/densityFunction.Eval(&x_i[i]);

        // for values below the defined range, the probability is zero
        if (x_trial < a)
            transitionProb = 0.0;

        if ((transitionProb >= 1.0 || transitionProb >= distribution(randomEngine))){

            x_i[i+1] = x_trial;
            acceptedSamples++;

        }

        else
            x_i[i+1] = x_i[i];


        if (acceptedSamples > 200){


            integralNumerator += (numeratorFunction.Eval(&x_i[i])/densityFunction.Eval(&x_i[i]))/samples;
            integralDenominator += (denominatorFunction.Eval(&x_i[i])/densityFunction.Eval(&x_i[i]))/samples;

            mean_value_numerator += numeratorFunction.Eval(&x_i[i])/samples;
            squared_mean_value_numerator += std::pow(numeratorFunction.Eval(&x_i[i]),2)/samples;


            mean_value_denominator += denominatorFunction.Eval(&x_i[i])/samples;
            squared_mean_value_denominator += std::pow(denominatorFunction.Eval(&x_i[i]),2)/samples;
        }

    }


    integralResult = integralNumerator/integralDenominator;


    standardDeviationNumerator = std::sqrt(squared_mean_value_numerator - std::pow(mean_value_numerator,2));
    standardDeviationDenominator = std::sqrt(squared_mean_value_denominator - std::pow(mean_value_denominator,2));


    this->stopClockAndPrintResults("Metropolis");


    delete [] x_i;

}




// Calculate the height of the rectangle, by taking the maximum value of the integrand
double NumericalMethods::getMaximumValueIntegrand(bool numerator){

    double x = a;
    double y_max, y_value;

    while(x < b){

        y_value = numerator ? numeratorFunction.Eval(&x) : denominatorFunction.Eval(&x);

        if (y_value > y_max)
            y_max = y_value;

        x += h;

    }

    return y_max;
}

// Calculate the height of the rectangle, by taking the maximum value of the integrand
double NumericalMethods::getXPositionForMaximumValue(){

    double x = a;
    double x_arg_max;
    double y_max = 0.0;

    while(x < b){

        if (densityFunction.Eval(&x) > y_max){
            y_max = densityFunction.Eval(&x);
            x_arg_max = x;

        }
        x += h;

    }

    return x_arg_max;
}

// Calculate standard deviation for method which slices the integral
void NumericalMethods::calculateStandardDeviation(){

    if (standardDeviationNumerator != 0.0 )
        return;

    double y_value_numerator, mean_value_numerator, mean_value_squared_numerator;
    double y_value_denominator, mean_value_denominator, mean_value_squared_denominator;
    double x_value = a;

    for (unsigned int i = 0; i < samples; i++){

        y_value_numerator = numeratorFunction.Eval(&x_value);
        y_value_denominator = denominatorFunction.Eval(&x_value);

        mean_value_numerator += y_value_numerator / samples;
        mean_value_squared_numerator += std::pow(y_value_numerator,2) / samples;

        mean_value_denominator += y_value_denominator / samples;
        mean_value_squared_denominator += std::pow(y_value_denominator,2) / samples;
        x_value += h;

    }

    standardDeviationNumerator = std::sqrt(mean_value_squared_numerator - std::pow(mean_value_numerator,2));
    standardDeviationDenominator = std::sqrt(mean_value_squared_denominator - std::pow(mean_value_denominator,2));

}

// Start the time and clear the integral value
void NumericalMethods::startClock() {

    timePassed = (double) clock();
    integralResult = 0.0;

}

// Stop clock and print the results
void NumericalMethods::stopClockAndPrintResults(std::string methodName) {

    timePassed = ((double) clock() - timePassed) * 1000.0 / CLOCKS_PER_SEC;

    this->calculateStandardDeviation();

    std::cout << methodName << " method:" << std::endl;
    std::cout << "Integral value: \t\t" << integralResult << std::endl;
    std::cout << "Error: \t\t\t\t\t" << std::abs(integralResult - analyticalSolution) << ", ";
    std::cout << std::abs(integralResult - analyticalSolution)*100/analyticalSolution << " %" << std::endl;
    std::cout << "Time: \t\t\t\t\t" << timePassed << " ms" << std::endl;

    // Note: The standard error is defined as \sigma / \sqrt{N}, with N being number of samples
    std::cout << "Standard error Num: \t" << standardDeviationNumerator/std::sqrt(samples) << ", ";
    std::cout << standardDeviationNumerator/std::sqrt(samples)*100 << " %" << std::endl;
    std::cout << "Standard error Den: \t" << standardDeviationDenominator/std::sqrt(samples) << ", ";
    std::cout << standardDeviationDenominator/std::sqrt(samples)*100 << " %" << std::endl << std::endl;

}

double NumericalMethods::getA() const {
    return a;
}

double NumericalMethods::getB() const {
    return b;
}

unsigned int NumericalMethods::getSamples() const {
    return samples;
}

double NumericalMethods::getAnalyticalSolution() const {
    return analyticalSolution;
}

