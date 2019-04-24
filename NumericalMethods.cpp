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

    timePassed = 0.0;
    integralResult = 0.0;

    stepSize = (b-a)/samples;

    // Define functions, parse from strings
    mainFunction.Parse(*stringFunctions[0], "x");
    weightFunction.Parse(*stringFunctions[1], "x");
    densityFunction.Parse(*stringFunctions[2], "x");

}


// Integral method using Simpson's method (integral slicing)
void NumericalMethods::simpson() {

    double numeratorResult = 0.0, denominatorResult = 0.0, x = a + stepSize, four_or_two;

    this->startClock();

    while (x < b) {

        // first line is either 4 or 2
        four_or_two = 2*((int)((x-a)/stepSize)%2+1);
        numeratorResult +=  four_or_two * mainFunction.Eval(&x) * weightFunction.Eval(&x) * 1/3.0 * stepSize;
        denominatorResult += four_or_two * weightFunction.Eval(&x) * 1/3.0 * stepSize;
        x += stepSize;

    }

    // add the boundary values of the Simpson's method
    numeratorResult += (mainFunction.Eval(&a)*weightFunction.Eval(&a) + mainFunction.Eval(&b)*weightFunction.Eval(&b)) * 1/3.0 * stepSize;
    denominatorResult += (weightFunction.Eval(&a) + weightFunction.Eval(&b)) * 1/3.0 * stepSize;

    integralResult = numeratorResult/denominatorResult;

    this->stopClockAndPrintResults("Simpson's");

}


// Integral calculation using Monte Carlo Mean Sample Value method
void NumericalMethods::simpleMonteCarlo(){

    auto * x_i = new double[samples];
    double integralNumerator = 0.0, integralDenominator = 0.0;

    this->startClock();

    // Create distribution for the range [0,1]
    std::random_device randomDevice;
    std::mt19937 randomEngine(randomDevice());
    std::uniform_real_distribution<> distribution(0.0, 1.0);

    for (int i = 0; i < samples; i++){

        // x randomized according to the PDF p(x) = 0.5exp(-x/2) from 0 to infinity -> x = -2ln(r)
        x_i[i] = -2.0*log(distribution(randomEngine));

        integralNumerator += (mainFunction.Eval(&x_i[i])*weightFunction.Eval(&x_i[i])/densityFunction.Eval(&x_i[i]));
        integralDenominator += (weightFunction.Eval(&x_i[i])/densityFunction.Eval(&x_i[i]));

    }

    integralResult = (integralNumerator/integralDenominator);

    this->stopClockAndPrintResults("Monte Carlo mean sample", this->calculateStandardError(x_i));

    delete [] x_i;
}

// Generate samples using Metropolis algorithm
void NumericalMethods::metropolis() {

    // Delta = size of random walk
    double delta = 10.0, integralNumerator = 0.0, integralDenominator = 0.0, x_trial, transitionProbability;
    auto * x_i = new double[samples];
    unsigned int acceptedSamples = 0;

    // Start the chain at maximum argument value
    x_i[0] = this->getXPositionForMaximumValue();

    // Create distribution for the range [0,1]
    std::random_device randomDevice;
    std::mt19937 randomEngine(randomDevice());
    std::uniform_real_distribution<> distribution(0.0, 1.0);

    this->startClock();

    for (unsigned int i = 0; i < samples; i++){

        x_trial = x_i[i] + (2.0 * distribution(randomEngine) - 1.0) * delta;
        transitionProbability = densityFunction.Eval(&x_trial)/densityFunction.Eval(&x_i[i]);

        // for values below the defined range, the probability is zero
        if (x_trial < a)
            transitionProbability = 0.0;

        if ((transitionProbability >= 1.0 || transitionProbability >= distribution(randomEngine))){

            x_i[i+1] = x_trial;
            acceptedSamples++;

        }
        else
            x_i[i+1] = x_i[i];


        if (acceptedSamples > 200){

            integralNumerator += (mainFunction.Eval(&x_i[i]) * weightFunction.Eval(&x_i[i])/densityFunction.Eval(&x_i[i]));
            integralDenominator += (weightFunction.Eval(&x_i[i])/densityFunction.Eval(&x_i[i]));

        }

    }

    integralResult = integralNumerator/integralDenominator;



    this->stopClockAndPrintResults("Metropolis", this->calculateStandardError(x_i));

    delete [] x_i;

}


// Calculate the height of the rectangle, by taking the maximum value of the integrand
double NumericalMethods::getXPositionForMaximumValue(){

    double x_value = a, y_max = 0.0, x_arg_max;

    while(x_value < b){

        if (densityFunction.Eval(&x_value) > y_max){
            y_max = densityFunction.Eval(&x_value);
            x_arg_max = x_value;

        }
        x_value += stepSize;

    }

    return x_arg_max;

}

// Calculate standard deviation for method which slices the integral
double NumericalMethods::calculateStandardError(double x_i[]){

    double y_value, x_value = a, mean_value = 0.0, mean_value_squared = 0.0;

    for (unsigned int i = 0; i < samples; i++){

        y_value = mainFunction.Eval(&x_i[i]);

        mean_value += y_value / samples;
        mean_value_squared += std::pow(y_value,2) / samples;
        x_value += stepSize;

    }

    return std::sqrt(mean_value_squared - std::pow(mean_value,2))/std::sqrt(samples);

}

// Start the time and clear the integral value
void NumericalMethods::startClock() {

    timePassed = (double) clock();
    integralResult = 0.0;

}

// Stop clock and print the results
void NumericalMethods::stopClockAndPrintResults(std::string methodName, double standardError) {

    timePassed = ((double) clock() - timePassed) * 1000.0 / CLOCKS_PER_SEC;

    std::cout << methodName << " method:" << std::endl;
    std::cout << "Integral value: \t\t" << integralResult << std::endl;
    std::cout << "Error: \t\t\t\t\t" << std::abs(integralResult - analyticalSolution) << ", ";
    std::cout << std::abs(integralResult - analyticalSolution)*100/analyticalSolution << " %" << std::endl;
    std::cout << "Time: \t\t\t\t\t" << timePassed << " ms" << std::endl;

    if (standardError != 0.0)
        std::cout << "Standard Error: \t\t" << standardError << ", " << standardError*100 << " %" << std::endl;

    std::cout << std::endl << std::endl;

}
