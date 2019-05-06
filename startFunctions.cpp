//
// Created by Antek Szadaj on 2019-03-13.
//

#include "startFunctions.h"
#include <string>
#include <sstream>

namespace plt = matplotlibcpp;


void obtainIntegralValue(double long * analyticalSolution, bool simpleIntegral){

    // Default sample = 1M, delta = 1.0
    NumericalMethods calculations = NumericalMethods(analyticalSolution, simpleIntegral);

    printIntegralInformation(&calculations, simpleIntegral);


    calculateIntegrals(&calculations);

    analyzeErrors(&calculations);

    analyzeCPUTimes(&calculations);

    analyzeCorrelationTime(&calculations);
    

}

void calculateIntegrals(NumericalMethods * calculations){

    std::cout << "Start usual integral calculation." << std::endl << std::endl;

    calculations->setSamples(1e6);

    calculations->delta = 1;

    calculations->simpson();
    calculations->simpleMonteCarlo();
    calculations->metropolis();

    std::cout << std::endl << "Done performing integral calculation." << std::endl << std::endl << std::endl;;

}


// Analyze the MC methods and calculate the RMS error
void analyzeErrors(NumericalMethods * calculations) {

    // the limit of increasing the number of samples, 10^{tenfold}
    int tenfold = 6;

    // Perform the RMS calculations per 10^N
    int numberOfIterationsPerSampleNumber = 100;

    calculations->delta = 10;

    double squared_deviation_error_metropolis, standard_error_average_metropolis;
    double squared_deviation_error_simple_monte_carlo, standard_error_average_simple_monte_carlo;
    int rejectedResultsMetropolis, rejectedResultsMonteCarlo;

    // 4 is rows here, number of different values
    std::vector<std::vector<double>> errors(4);
    for (int i = 0; i < 4; i++)
        errors[i].resize(tenfold);

    std::vector<double> x_ranges(tenfold);

    calculations->printMessage = false;

    std::cout << "Start error calculation." << std::endl << std::endl;

    std::cout << "Maximum walking step, aka delta: " << calculations->delta << std::endl;
    std::cout << "Samples: " << calculations->getSamples() << std::endl;
    std::cout << "Number of iterations per 10^N: " << numberOfIterationsPerSampleNumber << std::endl << std::endl;

    for (int j = 0; j <= tenfold-1; j++) {

        calculations->setSamples((int)std::pow(10, j+1));

        // reset values for each loop
        squared_deviation_error_metropolis = squared_deviation_error_simple_monte_carlo =
                standard_error_average_metropolis = standard_error_average_simple_monte_carlo = 0.0;

        rejectedResultsMetropolis = rejectedResultsMonteCarlo = 0;

        for (int i = 0; i < numberOfIterationsPerSampleNumber; i++) {

            calculations->metropolis();

            if (calculations->integralResult != 0.0) {

                squared_deviation_error_metropolis +=
                        std::pow((calculations->integralResult - calculations->analyticalSolution), 2);
                standard_error_average_metropolis +=
                        calculations->standardError / numberOfIterationsPerSampleNumber;
            }
            else
                rejectedResultsMetropolis++;

            calculations->simpleMonteCarlo();

            if (calculations->integralResult != 0.0) {
                squared_deviation_error_simple_monte_carlo +=
                        std::pow((calculations->integralResult - calculations->analyticalSolution), 2);
                standard_error_average_simple_monte_carlo +=
                        calculations->standardError / numberOfIterationsPerSampleNumber;

            }
            else
                rejectedResultsMonteCarlo++;

        }

        x_ranges[j] = calculations->getSamples();

        errors[0][j] = std::sqrt(squared_deviation_error_metropolis /
                (numberOfIterationsPerSampleNumber - rejectedResultsMetropolis));
        errors[1][j] = std::sqrt(squared_deviation_error_simple_monte_carlo /
                (numberOfIterationsPerSampleNumber - rejectedResultsMonteCarlo));
        errors[2][j] = standard_error_average_metropolis;
        errors[3][j] = standard_error_average_simple_monte_carlo;

        std::cout << "Samples: " << calculations->getSamples() << std::endl;
        std::cout << "Rejected Metropolis: " << rejectedResultsMetropolis << "/";
        std::cout << numberOfIterationsPerSampleNumber << std::endl;
        std::cout << "Rejected Monte Carlo: " << rejectedResultsMonteCarlo << "/";
        std::cout << numberOfIterationsPerSampleNumber << std::endl;
        std::cout << "RMS error Metropolis: " << errors[0][j] << std::endl;
        std::cout << "RMS error Simple Monte Carlo: " << errors[1][j] << std::endl;
        std::cout << "Average standard error Metropolis: " << errors[2][j] << std::endl;
        std::cout << "Average standard error Monte Carlo: " << errors[3][j] << std::endl << std::endl;

    }

    exportErrorPlot(&x_ranges, &errors, numberOfIterationsPerSampleNumber,
            calculations->delta, calculations->simpleIntegral);

    std::cout << std::endl << "Done performing error analysis." << std::endl << std::endl << std::endl;

}


void analyzeCPUTimes(NumericalMethods * calculations){

    std::cout << "Start CPU time analysis." << std::endl << std::endl;

    calculations->setSamples(1e6);
    calculations->delta = 2;

    if (calculations->simpleIntegral)
        calculations->errorLevel =1.0e-3;
    else
        calculations->errorLevel =1.0e-2;


    calculations->CPUTimeAnalysis = true;
    calculations->printMessage = true;

    calculations->simpson();

    std::cout << "Simpson's method: " << calculations->integralTime << " ms" << std::endl;
    std::cout << "Number of samples: " << calculations->sampleLevel << std::endl;
    std::cout << "Error level: " << calculations->errorLevel << std::endl << std::endl;

    calculations->simpleMonteCarlo();

    std::cout << "Simple Monte Carlo method: " << calculations->integralTime << " ms" << std::endl;
    std::cout << "Number of samples: " << calculations->sampleLevel << std::endl;
    std::cout << "Error level: " << calculations->errorLevel << std::endl << std::endl;

    calculations->metropolis();

    std::cout << std::endl << "Metropolis method: " << calculations->integralTime << " ms" << std::endl;
    std::cout << "Number of samples: " << calculations->sampleLevel << std::endl;
    std::cout << "Error level: " << calculations->errorLevel << std::endl << std::endl;

    std::cout << "Done with CPU time analysis." << std::endl << std::endl << std::endl;


}

void analyzeCorrelationTime(NumericalMethods * calculations){

    std::cout << "Start correlation time analysis." << std::endl << std::endl;

    calculations->setSamples((int)1e5);
    for (int i = 0; i < 3; i++){

        calculations->delta = 5*std::pow(10,i-1);
        calculations->metropolis();

    }

    calculations->setSamples((int)1e6);
    for (int i = 0; i < 3; i++){

        calculations->delta = 5*std::pow(10,i-1);
        calculations->metropolis();

    }

    std::cout << "Done with correlation time analysis." << std::endl << std::endl << std::endl;
}

// uses MATPLOTLIB to show the results.
void exportErrorPlot(std::vector<double> * x_ranges, std::vector<std::vector <double>> * errors,
        int numberOfIterationsPerSampleNumber, float delta, bool simpleIntegral){

    std::string exportDestination;

    if (simpleIntegral)
        exportDestination = "/Users/aszadaj/Desktop/SI2530 Computational Physics/Project/"
                                    "RMS_ASE_METROPOLIS_MC_FRAC_INTEGRAL.pdf";
    else
        exportDestination = "/Users/aszadaj/Desktop/SI2530 Computational Physics/Project/"
                            "RMS_ASE_METROPOLIS_MC_OSC_INTEGRAL.pdf";

    // this decreases the precision of digits for variable delta
    std::ostringstream out;
    out.precision(1);
    out << std::fixed << delta;
    out.str();

    // plot the figure with defined attributes
    plt::figure_size(1200, 780);
    plt::named_loglog("RMS Metropolis", *x_ranges, (*errors)[0], "b--o");
    plt::named_loglog("RMS MC", *x_ranges, (*errors)[1], "r--o");
    plt::named_loglog("Avg Standard Error Metropolis", *x_ranges, (*errors)[2], "g--d");
    plt::named_loglog("Avg Standard Error MC", *x_ranges, (*errors)[3], "m--d");
    std::string plot_title = "RMS and Average Standard Error for MC Methods, repeated "
            + std::to_string(numberOfIterationsPerSampleNumber) + " times for each N, $\\delta$ = " + out.str() + ".";
    plt::title(plot_title);
    plt::grid(true);
    plt::xlim(1.0, 1.0e7);
    plt::legend();

    plt::xlabel("Number of Samples [N]");
    plt::ylabel("Error (RMS or Average Standard Error)");

    plt::save(exportDestination);

    plt::close();

}


void printIntegralInformation(NumericalMethods * calculations, bool simpleIntegral){


    if (simpleIntegral){

        std::cout << std::endl << "Integral: int(x*exp(-x^2),";
        std::cout << calculations->getLowerLimit() << "," << calculations->getHigherLimit() << ",x)" << std::endl;
        std::cout << "Density function: 2*exp(-x^2)/sqrt(pi)"  << std::endl;
        std::cout << "Analytical solution: " << calculations->analyticalSolution;
    }
    else{

        std::cout << std::endl << "Integral: int(sin(1/x)^2,";
        std::cout << calculations->getLowerLimit() << "," << calculations->getHigherLimit() << ",x)" << std::endl;
        std::cout << "Density function: multiple"<< std::endl;
        std::cout << "Analytical solution: " << calculations->analyticalSolution;

    }


    std::cout << std::endl << std::endl << std::endl;

}
