
#include "NumericalMethods.h"

namespace plt = matplotlibcpp;

// Constructor
NumericalMethods::NumericalMethods(bool simpleIntegral) {

    integralResult = integralTime = standardError = errorFromRealResult = 0.0;
    delta = 1.0;
    samples = 1000000;
    errorLevel = 1.0e-4;
    rejectedSamples = 0;

    lowerLimit = 0.0001;
    higherLimit = 1000;

    printMessage = true;
    CPUTimeAnalysis = false;
    correlationTimeAnalysis = false;

    this->simpleIntegral = simpleIntegral;

    if (simpleIntegral)
        analyticalSolution = 0.5;
    else
        analyticalSolution = M_PI/2;

    stepSize = (higherLimit-lowerLimit)/samples;

}


// Integral method using Simpson's method
void NumericalMethods::simpson() {

    auto * x_i = new double[samples];
    float constant = 1/3.0 * stepSize;

    resetValues();
    startClock();

    // add the boundary values of the Simpson's method
    integralResult += (getMainFunction(&lowerLimit) + getMainFunction(&higherLimit)) * constant;

    for (unsigned int i = 1; i < samples; i++){

        x_i[i] = lowerLimit + i * stepSize;

        integralResult += 2*(i%2+1) * getMainFunction(&x_i[i]) * constant;

        // Check if the value reaches error level
        if (CPUTimeAnalysis && (std::abs(integralResult - analyticalSolution) < errorLevel))
            break;
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
    // Start engines to generate r = [0,1] and \theta = [0, pi/2]
    std::random_device thetaDevice, RDevice;
    std::mt19937 thetaEngine(thetaDevice()), REngine(RDevice());
    std::uniform_real_distribution<> theta_dist(0.0, M_PI/2.0), r_dist(0.0, 1.0);

    for (unsigned int i = 0; i < samples; i++){

        // x randomized based on \rho and \theta or the multiple choice pdf
        if (simpleIntegral)
            x_i[i] = std::sqrt(-std::log(r_dist(REngine))) *
                    std::cos(theta_dist(thetaEngine));

        else{

            nonUniformDistribution = r_dist(REngine);
            randomizedRValue = r_dist(REngine);
            randomizedNValue = r_dist(REngine);

            x_i[i] = getSampledPDFValue(&nonUniformDistribution, &randomizedRValue, &randomizedNValue);

        }

        integralResult += getMainFunction(&x_i[i]) / getPDF(&x_i[i]) / samples;

        // Check if the value reaches error level
        if (CPUTimeAnalysis && (std::abs(integralResult - analyticalSolution) < errorLevel))
            break;
    }

    stopClock();
    getStandardError(x_i);

    printResults("Simple Monte Carlo");

    delete [] x_i;

}

// Generate samples using Metropolis algorithm (The rejected samples are included in the integral)
void NumericalMethods::metropolis() {

    double x_trial, transProb;
    auto * x_i = new double[samples];

    resetValues();

    // Start the chain at maximum argument value
    x_i[0] = getRandomWalkStartValue();

    startClock();

    // Create distribution for the range [0,1]
    std::random_device device;
    std::mt19937 engine(device());
    std::uniform_real_distribution<> r_distr(0.0, 1.0);

    // Create random walk
    for (unsigned int i = 0; i < samples; i++){

        x_trial = x_i[i] + (2.0 * r_distr(engine) - 1.0) * delta;
        transProb = getPDF(&x_trial) / getPDF(&x_i[i]);

        // Imposed criterion on p(x) = 0 for x < a and x > b.
        if (x_trial < lowerLimit || x_trial > higherLimit)
            transProb = 0.0;

        if (transProb > 1.0 || transProb > r_distr(engine))
            x_i[i+1] = x_trial;

        else {
            x_i[i+1] = x_i[i];
            rejectedSamples++;
        }

        integralResult += getMainFunction(&x_i[i]) / getPDF(&x_i[i]) / samples;

        // Check if the value reaches a set error level
        if (CPUTimeAnalysis && (std::abs(integralResult - analyticalSolution) < errorLevel))
            break;
    }

    stopClock();
    getStandardError(x_i);

    if (correlationTimeAnalysis)
        exportCorrelationTimePlot(x_i);

    printResults("Metropolis");

    delete [] x_i;

}


// Given an array with size N, create a chain of x-values according to q(x)
void NumericalMethods::createRandomWalk(double x_i[]){

    double x_trial, transitionProbability;

    // Start randomization engine, to generate r on [0,1]
    std::random_device randomDevice;
    std::mt19937 randomEngine(randomDevice());
    std::uniform_real_distribution<> distribution(0.0, 1.0);

    // Obtain maximum value to start the random walk
    x_i[0] = getRandomWalkStartValue();

    for (unsigned int i = 0; i < samples; i++){

        x_trial = x_i[i] + (2.0 * distribution(randomEngine) - 1.0) * delta;
        transitionProbability = getPDF(&x_trial) / getPDF(&x_i[i]);

        // Imposed criterion on q(x) = 0 for x < a and x > b.
        if (x_trial < lowerLimit || x_trial > higherLimit)
            transitionProbability = 0.0;

        if (transitionProbability > 1.0 ||
            transitionProbability > distribution(randomEngine))
            x_i[i+1] = x_trial;

        else
            x_i[i+1] = x_i[i];
    }
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


void NumericalMethods::exportCorrelationTimePlot(double x_i []){

    std::string exportDestination;
    std::vector<std::string> markers(3);
    markers[0] = "b--o";
    markers[1] = "r--o";
    markers[2] = "g--o";
    plt::figure_size(1200, 780);

    // this decreases the precision of digits for variable delta
    std::ostringstream out;

    int numberOfValues = 50;

    std::vector<float> correlationValues (numberOfValues);
    std::vector<float> x_values (numberOfValues);

    for (int i = 0; i < 3; i++){

        delta = std::pow(10,i-1);

//        if (i == 1)
//            numberOfValues = 20;
//        else
//            numberOfValues = 50;


        createRandomWalk(x_i);

        std::cout << std::endl << "delta: " << delta << std::endl;


        for (int i = 0; i < numberOfValues; i++){

            correlationValues[i] = getAutocorrelationValue(x_i, i);
            x_values[i] = i;

        }

        out.precision(1);
        out << std::fixed << delta;

        plt::title("Autocorrelation function values for various $\\delta$");
        plt::xlabel("Sample j");
        plt::ylabel("C(j)");
        plt::named_semilogy( "$\\delta$ = " + out.str(), x_values, correlationValues, markers[i]);
        out.str("");
        std::fill(correlationValues.begin(), correlationValues.end(), 0.0);
        std::fill(x_values.begin(), x_values.end(), 0.0);

    }

    exportDestination = "/Users/aszadaj/Desktop/SI2530 Computational Physics/Project/report/figures/"
                        "distribution_correlationTime_"+getFunctionName()+".pdf";

    plt::grid(true);
    plt::legend();
    plt::xlim(1.0, 80.0);

    plt::save(exportDestination);
    plt::close();

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
            functionValue = normalizationConstant * std::exp(2/M_PI - *x);

        else if (*x < 2/M_PI && *x > 1/M_PI)
            functionValue = normalizationConstant * (*x - 1/M_PI);

        else{

            waveNumberN = (int) (1.0/(M_PI**x));
            functionValue = normalizationConstant *
                    std::sin( waveNumberN * std::pow(M_PI,2) * (waveNumberN + 1)**x - waveNumberN * M_PI );
        }
    }

    return functionValue;
}


double NumericalMethods::getSampledPDFValue(double * randomizedPartOfPDF, double * randomizedRValue,
        double * randomizedNValue){

    double waveN, sampledXValue;
    double exponentialLimit = 2*std::pow(M_PI,2)/(4+M_PI+2*std::pow(M_PI,2));
    double linearLimit = exponentialLimit + M_PI/(4+M_PI+2*std::pow(M_PI,2));


    if (*randomizedPartOfPDF < exponentialLimit)
        sampledXValue = 2/M_PI - std::log(*randomizedRValue);

    else if(*randomizedPartOfPDF < linearLimit)
        sampledXValue = 1/M_PI + std::sqrt(*randomizedRValue);

    else{
        waveN = (int)(1.0/ *randomizedNValue);
        sampledXValue = (std::acos(1.0-2.0**randomizedRValue) + waveN * M_PI)/
                (waveN * std::pow(M_PI,2) * (waveN + 1));
    }
    return sampledXValue;
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

    integralResult = integralTime = errorFromRealResult = standardError =  0.0;
    rejectedSamples = 0;

}


unsigned int NumericalMethods::getSamples() const {
    return samples;
}


void NumericalMethods::setSamples(unsigned int samples) {

    NumericalMethods::samples = samples;
    NumericalMethods::stepSize = (higherLimit-lowerLimit)/samples;

}


double NumericalMethods::getLowerLimit() const {
    return lowerLimit;
}


double NumericalMethods::getHigherLimit() const {
    return higherLimit;
}


std::string NumericalMethods::getFunctionName(){


    if (simpleIntegral)
        return "simpleIntegral";
    else
        return "oscillatory";

}