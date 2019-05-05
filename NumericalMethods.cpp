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
NumericalMethods::NumericalMethods(const double long * analyticalSolution, bool simpleIntegral) {

    integralResult = integralTime = standardError = errorFromRealResult = correlationTime = 0.0;
    delta = 1.0;
    samples = sampleLevel = 1000000;
    errorLevel = 1.0e-4;
    rejectedSamples = 0;

    lowerLimit = 0.001;
    higherLimit = 100;

    printMessage = true;
    CPUTimeAnalysis = false;

    this->simpleIntegral = simpleIntegral;
    this->analyticalSolution = *analyticalSolution;

    stepSize = (higherLimit-lowerLimit)/samples;

}


// Integral method using Simpson's method
void NumericalMethods::simpson() {

    auto * x_i = new double[samples];
    double simpsonsConstant = 1/3.0 * stepSize;

    resetValues();
    startClock();

    // add the boundary values of the Simpson's method
    integralResult += (getMainFunction(&lowerLimit) + getMainFunction(&higherLimit)) * simpsonsConstant;

    for (unsigned int i = 1; i < samples-1; i++){

        x_i[i] = lowerLimit + i * stepSize;

        integralResult += 2*(i%2+1) * getMainFunction(&x_i[i]) * simpsonsConstant;

        // Check if the value reaches error level
        if (CPUTimeAnalysis && (std::abs(integralResult - analyticalSolution) < errorLevel)){
            sampleLevel = i;
            break;
        }
    }

    stopClock();
    printResults("Simpson's");

    delete [] x_i;

}


// Integral calculation using Simple Monte Carlo method
void NumericalMethods::simpleMonteCarlo(){

    auto * x_i = new double[samples];
    double nonUniformDistribution, randomizedRValue, randomizedNValue;

    resetValues();
    startClock();

    // Sample \rho on sqrt(-ln(r)) with r on [0,1] and \theta on [0, pi/2] for randomizing x
    std::random_device randomDevice_theta, randomDevice_R;
    std::mt19937 randomEngine_theta(randomDevice_theta()), randomEngine_R(randomDevice_R());
    std::uniform_real_distribution<> distribution_theta(0.0, M_PI/2.0), distribution_R(0.0, 1.0);

    for (unsigned int i = 0; i < samples; i++){

        // x randomized based on \rho and \theta or p(x) = 0.5*x*exp(-x^2)
        if (simpleIntegral)
            x_i[i] = std::sqrt(-std::log(distribution_R(randomEngine_R))) *
                    std::cos(distribution_theta(randomEngine_theta));

        else{

            nonUniformDistribution = distribution_R(randomEngine_R);
            randomizedRValue = distribution_R(randomEngine_R);
            randomizedNValue = distribution_R(randomEngine_R);

            x_i[i] = getSampledPDFValue(&nonUniformDistribution, &randomizedRValue, &randomizedNValue);

        }

        integralResult += getMainFunction(&x_i[i]) / getPDF(&x_i[i]) / samples;

        // Check if the value reaches error level
        if (CPUTimeAnalysis && (std::abs(integralResult - analyticalSolution) < errorLevel)){
            sampleLevel = i;
            break;
        }
    }

    stopClock();
    getStandardError(x_i);

    printResults("Simple Monte Carlo");
    //exportRandomizedSamples(x_i, "monte_carlo");

    delete [] x_i;

}

// Generate samples using Metropolis algorithm (warning, the rejected samples are included in the integral)
void NumericalMethods::metropolis() {

    double transitionProbability, x_trial;
    auto * x_i = new double[samples];

    resetValues();

    // Start the chain at maximum argument value
    x_i[0] = getRandomWalkStartValue();

    startClock();

    // Create distribution for the range [0,1]
    std::random_device randomDevice;
    std::mt19937 randomEngine(randomDevice());
    std::uniform_real_distribution<> distribution(0.0, 1.0);

    // Create the random walk
    for (unsigned int i = 0; i < samples; i++){

        x_trial = x_i[i] + (2.0 * distribution(randomEngine) - 1.0) * delta;
        transitionProbability = getPDF(&x_trial) / getPDF(&x_i[i]);

        // Imposed criterion on p(x) = 0 for x < a and x > b.
        if (x_trial < lowerLimit || x_trial > higherLimit)
            transitionProbability = 0.0;

        if (transitionProbability > 1.0 || transitionProbability > distribution(randomEngine))
            x_i[i+1] = x_trial;

        else {
            x_i[i+1] = x_i[i];
            rejectedSamples++;
        }

        integralResult += getMainFunction(&x_i[i]) / getPDF(&x_i[i]) / samples;

        // Check if the value reaches a set error level
        if (CPUTimeAnalysis && (std::abs(integralResult - analyticalSolution) < errorLevel)){
            sampleLevel = i;
            break;
        }
    }

    stopClock();

    getStandardError(x_i);

    // the result is in us (microseconds)
    correlationTime = integralTime * (float) determineCorrelatedStep(x_i)/samples * 1000;

    printResults("Metropolis");
    //exportRandomizedSamples(x_i, "metropolis");

    delete [] x_i;

}


double NumericalMethods::getRandomWalkStartValue(){

    double x_value = lowerLimit, x_arg_max = lowerLimit, y_value, y_max = 0.0;

    while (x_value < higherLimit){

        y_value = getPDF(&x_value);

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

        y_value = getMainFunction(&x_i[i])/getPDF(&x_i[i]);

        mean_value += y_value / samples;
        mean_value_squared += std::pow(y_value,2) / samples;
        x_value += stepSize;

    }

    standardError = std::sqrt(mean_value_squared - std::pow(mean_value,2))/std::sqrt(samples);

}



unsigned int NumericalMethods::determineCorrelatedStep(double x_i[]){

    double autocorrelationValue;
    float correlationLimit = 0.0;
    unsigned int correlatedStep = samples;

    for (unsigned int j = 0; j < samples; j++) {

        autocorrelationValue = getAutocorrelationValue(x_i, j);

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


double NumericalMethods::getMainFunction(double * x){

    if (simpleIntegral)
        return *x * std::exp(-std::pow(*x,2));
    else
        return std::pow(std::sin(1.0 / *x),2);


}


double NumericalMethods::getPDF(double * x){

    // wave number on the oscillatory part, counted down towards x = 0.
    int waveNumberN;
    double functionValue;
    double normalizationConstant = 2*std::pow(M_PI,2) / (4 + M_PI + 2*std::pow(M_PI,2));

    if (simpleIntegral)
        functionValue = 2.0*std::exp(-std::pow(*x,2))/std::sqrt(M_PI);

    else{

        if (*x >= 2/M_PI)
            functionValue = normalizationConstant*std::exp(2/M_PI - *x);

        else if (*x < 2/M_PI && *x > 1/M_PI)
            functionValue = normalizationConstant*(*x - 1/M_PI);

        else{

            waveNumberN = (int) (1.0/(M_PI**x));
            functionValue = normalizationConstant *
                    std::sin( waveNumberN * std::pow(M_PI,2) * (waveNumberN + 1)**x - waveNumberN * M_PI );
        }
    }

    return functionValue;
}


double NumericalMethods::getSampledPDFValue(double * nonUniformDistribution, double * randomizedRValue,
        double * randomizedNValue ){

    double waveNumberN, sampledOnRangeXValue, oscillatoryPartValue;
    double exponentialLimit = 2*std::pow(M_PI,2)/(4+M_PI+2*std::pow(M_PI,2));
    double linearLimit = exponentialLimit + M_PI/(4+M_PI+2*std::pow(M_PI,2));


    if (*nonUniformDistribution < exponentialLimit){
        sampledOnRangeXValue = 2/M_PI - std::log(*randomizedRValue);
    }


    else if(*nonUniformDistribution < linearLimit && *nonUniformDistribution > exponentialLimit){
        sampledOnRangeXValue = 1/M_PI + std::sqrt(*randomizedRValue);
    }


    else{

        waveNumberN = (int)(1.0/ *randomizedNValue);
        sampledOnRangeXValue = (std::acos(1.0-2.0**randomizedRValue) + waveNumberN * M_PI)/
                (waveNumberN * std::pow(M_PI,2) * (waveNumberN + 1));
    }

    return sampledOnRangeXValue;
}




void NumericalMethods::exportRandomizedSamples(double x_i[], std::string function){


    std::string exportDestination = "/Users/aszadaj/Desktop/SI2530 Computational Physics/Project/"
                                    "distribution_"+function+".pdf";


    std::vector<double> x(samples);
    for (int i = 0; i < samples; i++)
        x[i] = x_i[i];

    plt::hist(x, 900);
    //plt::xlim(0, 5);

    // plot the figure with defined attributes
    //plt::figure_size(1200, 780);

    plt::save(exportDestination);

    plt::close();
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
    std::cout << "Standard Error: \t\t" << standardError << ", ";
    std::cout << standardError*100 << " %" << std::endl;
    std::cout << "Correlation time: \t\t" << correlationTime << " us" << std::endl;
    std::cout << "Rejected samples:\t\t" << (float)rejectedSamples/samples*100<< " %" << std::endl;

    std::cout << std::endl << std::endl;

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


void NumericalMethods::resetValues(){

    correlationTime  = integralResult = integralTime = errorFromRealResult = standardError =  0.0;
    rejectedSamples = 0;
    sampleLevel = samples;

}


unsigned int NumericalMethods::getSamples() const {
    return samples;
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
