//
// Created by Antek Szadaj on 2019-02-19.
//

#include "NumericalMethods.h"

// Constructor
NumericalMethods::NumericalMethods(double * parameters [], std::string * stringFunctions []) {

    timePassed = 0.0;
    integralResult = 0.0;
    standardError = 0.0;
    errorFromRealResult = 0.0;
    correlationTime = 0.0;
    printMessage = true;

    // introduce 0.001 % level to determine how fast the method reaches this level
    CPUTimeAnalysis = false;
    errorLevel = 1.0e-5;

    // sample level is the level at which the error level is reached. Default the maximum.
    samples = sampleLevel = 1000000;

    lowerLimit = 0.0;
    higherLimit = 1.0;

    analyticalSolution = *parameters[0];

    stepSize = (higherLimit-lowerLimit)/samples;
    delta = *parameters[1];
    randomWalkTime = 0.0;

    // Define functions, parse from strings
    integrandFunction.Parse(*stringFunctions[0], "x");
    densityFunction.Parse(*stringFunctions[1], "x");

}


// Integral method using Simpson's method
void NumericalMethods::simpson() {

    auto * x_i = new double [samples];
    double four_or_two, simpsonsConstant = 1/3.0 * stepSize;

    this->resetValues();
    this->startClock();

    // add the boundary values of the Simpson's method
    integralResult += (integrandFunction.Eval(&lowerLimit) + integrandFunction.Eval(&higherLimit)) * simpsonsConstant;

    for (unsigned int i = 1; i < samples-1; i++){

        x_i[i] = lowerLimit + i * stepSize;
        four_or_two = 2*(i%2+1);
        integralResult += four_or_two * integrandFunction.Eval(&x_i[i]) * simpsonsConstant;

        if (CPUTimeAnalysis && (std::abs(integralResult - analyticalSolution) < errorLevel)){
            sampleLevel = i;
            break;
        }

    }

    this->stopClock();
    this->getStandardError(x_i, false);
    this->printResults("Simpson's");

    delete [] x_i;

}


// Integral calculation using Monte Carlo Mean Sample Value method
void NumericalMethods::simpleMonteCarlo(){

    auto * x_i = new double[samples];
    this->resetValues();
    this->startClock();

    // Create distribution for the range [0,1] and [0, pi/2]
    // for creating distributions for \rho and \theta
    std::random_device randomDevice_theta, randomDevice_R;
    std::mt19937 randomEngine_theta(randomDevice_theta()), randomEngine_R(randomDevice_R());
    std::uniform_real_distribution<> distribution_theta(0.0, M_PI/2.0), distribution_R(0.0, 1.0);

    /*

     warning!

     since p(x) = 2/sqrt(pi) * exp(-x^2) cannot be solved for x analytically, impose the box muller method

     generate two gaussian distributions, such that

     p(\rho,\theta) d\rho d\theta = 2/sqrt(pi) exp(-\rho) d\rho d\theta

     is obtained. Therefore, the distributions are as follows

     x = \sqrt(\rho) \cos \theta,
     y = \sqrt(\rho) \sin \theta,

     where \rho is distributed according to \exp(-\rho), which yields \rho = -\ln(r),
     with r uniformly distributed on [0,1]

     \theta is uniformly distributed on [0,\pi/2]

     this can be seen directly from

     p(\rho,\theta) d\rho d\theta = 2/sqrt(pi) d\theta exp(-\rho) d\rho
                                    _________________  ________________
                                          = \theta           = \rho

    */

    for (unsigned int i = 0; i < samples; i++){

        // x randomized based on \rho and \theta
//        x_i[i] = std::sqrt(-log(distribution_R(randomEngine_R))) * cos(distribution_theta(randomEngine_theta));

        x_i[i] = -log(1.0 - distribution_R(randomEngine_R)*(exp(1)-1.0)/exp(1));

        integralResult += (integrandFunction.Eval(&x_i[i]) / densityFunction.Eval(&x_i[i]))/samples;

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

    unsigned int correlatedStep;
    auto * x_i = new double[samples];

    this->resetValues();
    this->createRandomWalk(x_i);

    correlatedStep = this->determineCorrelatedStep(x_i);

    // Exclude the time it takes to find the correlated step, and adjust for the time it takes to do the random walk
    this->startClock(randomWalkTime);

    for (unsigned int i = correlatedStep; i < samples; i++){

        integralResult += integrandFunction.Eval(&x_i[i]) / densityFunction.Eval(&x_i[i]) / samples;

        if (CPUTimeAnalysis && (std::abs(integralResult - analyticalSolution) < errorLevel)){
            sampleLevel = i;
            break;
        }
    }

    this->stopClock();
    this->getStandardError(x_i);

    // the result is in us (microseconds)
    correlationTime = randomWalkTime * (float) correlatedStep/samples * 1000;

    this->printResults("Metropolis");

    delete [] x_i;

}


// Calculate the height of the rectangle, by taking the maximum value of the integrand
double NumericalMethods::getArgMax(){

    double x_value = lowerLimit, x_arg_max = lowerLimit, y_value, y_max = 0.0;

    while (x_value < higherLimit){

        y_value = densityFunction.Eval(&x_value);

        if (y_value > y_max){
            y_max = y_value;
            x_arg_max = x_value;
        }

        x_value += stepSize;

    }

    return x_arg_max;

}

// Calculate standard error based on the main function
void NumericalMethods::getStandardError(double x_i[], bool MC_METHODS){

    double y_value, x_value = lowerLimit, mean_value = 0.0, mean_value_squared = 0.0;

    for (unsigned int i = 0; i < samples; i++){

        y_value = MC_METHODS ? integrandFunction.Eval(&x_i[i])/densityFunction.Eval(&x_i[i]) : integrandFunction.Eval(&x_i[i]);

        mean_value += y_value / samples;
        mean_value_squared += std::pow(y_value,2) / samples;
        x_value += stepSize;

    }

    standardError = std::sqrt(mean_value_squared - std::pow(mean_value,2))/std::sqrt(samples);

}

// Start the time
void NumericalMethods::startClock(double additionalTime) {

    timePassed = (double) clock() - additionalTime * CLOCKS_PER_SEC/1000.0;

}

// Stop clock and print the results
void NumericalMethods::stopClock() {

    timePassed = ((double) clock() - timePassed) * 1000 / CLOCKS_PER_SEC;
    errorFromRealResult = std::abs(integralResult - analyticalSolution);

}


void NumericalMethods::createRandomWalk(double x_i[]){

    double transitionProbability, x_trial;

    // Start the chain at maximum argument value
    x_i[0] = this->getArgMax();

    auto startTime = (double) clock();

    // Create distribution for the range [0,1]
    std::random_device randomDevice;
    std::mt19937 randomEngine(randomDevice());
    std::uniform_real_distribution<> distribution(0.0, 1.0);

    for (unsigned int i = 0; i < samples; i++){

        x_trial = x_i[i] + (2 * distribution(randomEngine) - 1) * delta;
        transitionProbability = densityFunction.Eval(&x_trial) / densityFunction.Eval(&x_i[i]);

        // for values below the defined range, the probability is zero
        if (x_trial < lowerLimit)
            transitionProbability = 0.0;

        // Accept or reject walking point
        x_i[i + 1] =
                transitionProbability >= 1.0 || transitionProbability >= distribution(randomEngine) ? x_trial : x_i[i];

    }

    // result in ms
    randomWalkTime = ((double) clock() - startTime) * 1000 / CLOCKS_PER_SEC;


}

unsigned int NumericalMethods::determineCorrelatedStep(double x_i[]){

    double autocorrelationValue;
    float correlationLimit = 0.01;
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
    std::cout << "Time: \t\t\t\t\t" << timePassed << " ms" << std::endl;
    std::cout << "Standard Error: \t\t" << standardError << ", " << standardError*100 << " %" << std::endl;
    std::cout << "Correlation time: \t\t" << correlationTime << " us" << std::endl;
    std::cout << "Random walk time: \t\t" << randomWalkTime << " ms" << std::endl;

    std::cout << std::endl << std::endl;

}

void NumericalMethods::setSamples(unsigned int samples) {

    NumericalMethods::samples = samples;
    NumericalMethods::sampleLevel = samples;
    NumericalMethods::stepSize = (higherLimit-lowerLimit)/samples;

}

void NumericalMethods::resetValues(){

    correlationTime = randomWalkTime = integralResult = timePassed = errorFromRealResult = standardError =  0.0;

}

unsigned int NumericalMethods::getSamples() const {
    return samples;
}
