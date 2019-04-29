//
// Created by Antek Szadaj on 2019-02-19.
//

/*

 warning (regarding simpleMonteCarlo()-method)

 since p(x) = 2/sqrt(pi) * exp(-x^2) cannot be solved for x analytically, impose the box muller method
 generate two gaussian distributions, such that

    p(\rho,\theta) d\rho d\theta = 2/sqrt(pi) exp(-\rho) d\rho d\theta

 is obtained. Therefore, the distributions are as follows

     x = \sqrt(\rho) \cos \theta,
     y = \sqrt(\rho) \sin \theta,

 where \rho is distributed according to \exp(-\rho), which yields \rho = -\ln(r),
 with r uniformly distributed on [0,1] and \theta is uniformly distributed on [0,\pi/2].
 This can be seen directly from

    p(\rho,\theta) d\rho d\theta = 2/sqrt(pi) d\theta exp(-\rho) d\rho
                                   __________________ ________________
                                       = \theta            = \rho
*/


#include "NumericalMethods.h"

namespace plt = matplotlibcpp;

// Constructor
NumericalMethods::NumericalMethods(const double long * analyticalSolution, float * ranges[], std::string * stringFunctions []) {

    integralResult = integralTime = standardError = errorFromRealResult = randomWalkTime = correlationTime = 0.0;
    delta = 1.0;
    samples = sampleLevel = 1000000;
    errorLevel = 1.0e-4;
    rejectedSamples = 0;

    printMessage = true;
    CPUTimeAnalysis = false;

    this->analyticalSolution = *analyticalSolution;
    lowerLimit = *ranges[0];
    higherLimit = *ranges[1];

    // Define functions, parse from strings
    integrandFunction.Parse(*stringFunctions[0], "x");
    PDFFunction.Parse(*stringFunctions[1], "x");

    stepSize = (higherLimit-lowerLimit)/samples;

}


// Integral method using Simpson's method
void NumericalMethods::simpson() {

    double simpsonsConstant = 1/3.0 * stepSize;

    auto * x_i = new (std::nothrow) double[samples];

    if (x_i == nullptr){
        std::cout << "ERROR: Memory could not be allocated." << std::endl << std::endl;
        return;
    }

    this->resetValues();

    this->startClock();

    // add the boundary values of the Simpson's method
    integralResult += (integrandFunction.Eval(&lowerLimit) + integrandFunction.Eval(&higherLimit)) * simpsonsConstant;

    for (unsigned int i = 1; i < samples-1; i++){

        x_i[i] = lowerLimit + i * stepSize;

        integralResult += 2*(i%2+1) * integrandFunction.Eval(&x_i[i]) * simpsonsConstant;

        // Check if the value reaches error level
        if (CPUTimeAnalysis && (std::abs(integralResult - analyticalSolution) < errorLevel)){
            sampleLevel = i;
            break;
        }

    }

    this->stopClock();

    this->printResults("Simpson's");

    delete [] x_i;

}


// Integral calculation using Simple Monte Carlo method
void NumericalMethods::simpleMonteCarlo(){

    auto * x_i = new double[samples];

    this->resetValues();
    this->startClock();

    // Sample \rho on sqrt(-ln(r)) with r on [0,1] and \theta on [0, pi/2] for randomizing x
    std::random_device randomDevice_theta, randomDevice_R;
    std::mt19937 randomEngine_theta(randomDevice_theta()), randomEngine_R(randomDevice_R());
    std::uniform_real_distribution<> distribution_theta(0.0, M_PI/2.0), distribution_R(0.0, 1.0);

    for (unsigned int i = 0; i < samples; i++){

        // x randomized based on \rho and \theta
        x_i[i] = std::sqrt(-log(distribution_R(randomEngine_R))) * cos(distribution_theta(randomEngine_theta));

        integralResult += (integrandFunction.Eval(&x_i[i]) / PDFFunction.Eval(&x_i[i])) / samples;

        // Check if the value reaches error level
        if (CPUTimeAnalysis && (std::abs(integralResult - analyticalSolution) < errorLevel)){
            sampleLevel = i;
            break;
        }
    }

    this->stopClock();

    this->getStandardError(x_i);

    this->printResults("Simple Monte Carlo");

    delete [] x_i;

}

// Generate samples using Metropolis algorithm
void NumericalMethods::metropolis() {

    double transitionProbability, x_trial;
    auto * x_i = new double[samples];

    this->resetValues();

    // Start the chain at maximum argument value
    x_i[0] = this->getArgMax();

    this->startClock();

    // Create distribution for the range [0,1]
    std::random_device randomDevice;
    std::mt19937 randomEngine(randomDevice());
    std::uniform_real_distribution<> distribution(0.0, 1.0);

    // Create the random walk
    for (unsigned int i = 0; i < samples; i++){

        x_trial = x_i[i] + (2.0 * distribution(randomEngine) - 1.0) * delta;
        transitionProbability = PDFFunction.Eval(&x_trial) / PDFFunction.Eval(&x_i[i]);

        if (x_trial < lowerLimit || x_trial > higherLimit)
            transitionProbability = 0.0;

        if (transitionProbability > 1.0 || transitionProbability > distribution(randomEngine))
            x_i[i+1] = x_trial;

        else {
            x_i[i+1] = x_i[i];
            rejectedSamples++;
        }

        integralResult += integrandFunction.Eval(&x_i[i]) / PDFFunction.Eval(&x_i[i]) / samples;

        // Check if the value reaches a set error level
        if (CPUTimeAnalysis && (std::abs(integralResult - analyticalSolution) < errorLevel)){
            sampleLevel = i;
            break;
        }
    }

    this->stopClock();

    this->getStandardError(x_i);

    // the result is in us (microseconds)
    correlationTime = integralTime * (float) this->determineCorrelatedStep(x_i)/samples * 1000;

    this->printResults("Metropolis");

    delete [] x_i;

}


double NumericalMethods::getArgMax(){

    double x_value = lowerLimit, x_arg_max = lowerLimit, y_value, y_max = 0.0;

    while (x_value < higherLimit){

        y_value = PDFFunction.Eval(&x_value);

        if (y_value > y_max){
            y_max = y_value;
            x_arg_max = x_value;
        }

        x_value += stepSize;

    }

    return x_arg_max;

}


// Calculate standard error based on the main function
void NumericalMethods::getStandardError(double x_i[]){

    double y_value, x_value = lowerLimit, mean_value = 0.0, mean_value_squared = 0.0;

    for (unsigned int i = 0; i < samples; i++){

        y_value = integrandFunction.Eval(&x_i[i])/PDFFunction.Eval(&x_i[i]);

        mean_value += y_value / samples;
        mean_value_squared += std::pow(y_value,2) / samples;
        x_value += stepSize;

    }

    standardError = std::sqrt(mean_value_squared - std::pow(mean_value,2))/std::sqrt(samples);

}


// Start the time
void NumericalMethods::startClock(double additionalTime) {

    integralTime = (double) clock() - additionalTime * CLOCKS_PER_SEC/1000.0;

}


// Stop clock and print the results
void NumericalMethods::stopClock() {

    integralTime = ((double) clock() - integralTime) * 1000 / CLOCKS_PER_SEC;
    errorFromRealResult = std::abs(integralResult - analyticalSolution);

}


unsigned int NumericalMethods::determineCorrelatedStep(double x_i[]){

    double autocorrelationValue;
    float correlationLimit = 0.0;
    unsigned int correlatedStep = samples;

    for (unsigned int j = 0; j < samples; j++) {

        autocorrelationValue = this->getAutocorrelationValue(x_i, j);

        if (autocorrelationValue < correlationLimit) {
            correlatedStep = j;
            break;
        }
    }

    return correlatedStep;
}


double NumericalMethods::getAutocorrelationValue(double x_i [], int j){

    double mean_value = 0.0, mean_value_squared = 0.0, mean_value_convoluted = 0.0, autocorrelationValue;

    for (unsigned int k = 0; k < samples - j; k++){

        mean_value += x_i[k] / (samples - j);
        mean_value_squared += std::pow(x_i[k],2) / (samples - j);
        mean_value_convoluted += x_i[k] * x_i[k+j] / (samples - j);

    }

    autocorrelationValue = (mean_value_convoluted - std::pow(mean_value,2))/
            (mean_value_squared - std::pow(mean_value,2));

    return autocorrelationValue;

}


void NumericalMethods::printResults(std::string methodName){

    if (not printMessage)
        return;

    std::cout << methodName << " method:" << std::endl;
    std::cout << "Samples:\t\t\t\t" << samples << std::endl;
    std::cout << "Delta:\t\t\t\t\t" << delta << std::endl;
    std::cout << "Integral value: \t\t" << integralResult << std::endl;
    std::cout << "Error: \t\t\t\t\t" << errorFromRealResult << std::endl;
    std::cout << "Time: \t\t\t\t\t" << integralTime << " ms" << std::endl;
    std::cout << "Standard Error: \t\t" << standardError << ", " << standardError*100 << " %" << std::endl;
    std::cout << "Correlation time: \t\t" << correlationTime << " us" << std::endl;
    std::cout << "Random walk time: \t\t" << randomWalkTime << " ms" << std::endl;
    std::cout << "Rejected samples:\t\t" << (float)rejectedSamples/samples*100<< " %" << std::endl;

    std::cout << std::endl << std::endl;

}


void NumericalMethods::resetValues(){

    correlationTime = randomWalkTime = integralResult = integralTime = errorFromRealResult = standardError =  0.0;
    rejectedSamples = 0;
    sampleLevel = samples;

}


unsigned int NumericalMethods::getSamples() const {
    return samples;
}


void NumericalMethods::setLowerLimit(double lowerLimit) {
    NumericalMethods::lowerLimit = lowerLimit;
    this->setSamples(samples);
}


void NumericalMethods::setHigherLimit(double higherLimit) {
    NumericalMethods::higherLimit = higherLimit;
    this->setSamples(samples);
}


void NumericalMethods::setSamples(unsigned int samples) {

    NumericalMethods::samples = samples;
    NumericalMethods::sampleLevel = samples;
    NumericalMethods::stepSize = (higherLimit-lowerLimit)/samples;

}


double NumericalMethods::getLowerLimit() const {
    return lowerLimit;
}


double NumericalMethods::getHigherLimit() const {
    return higherLimit;
}
