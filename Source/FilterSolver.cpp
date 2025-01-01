/*
  ==============================================================================

    FilterSolver.cpp
    Created: 30 Dec 2024 7:47:50pm
    Author:  Barak

  ==============================================================================
*/

#include "FilterSolver.h"
#include <iostream>

FilterSolver::FilterSolver(std::vector<int> drawnPoints, double sampleRate) :
    m_drawnPoints(drawnPoints),
    m_sampleRate(sampleRate)
{
}
FilterSolver::~FilterSolver() {}

double FilterSolver::solve(const dlib::matrix<double, 3 * NUM_FILTERS, 1>& variables) {

    double difference = 0.f;

    updateCoefficientsBulk(variables);

    // Go through each drawn point, convert it to a frequency, get the filter response
    // and compare to expected response
    for (int i = 0; i < m_drawnPoints.size(); ++i) {
        double mag = 1.f;
        double pointToFreq = juce::mapToLog10(float(i / m_drawnPoints.size()), FREQ_RANGE.first, FREQ_RANGE.second);

        mag *= m_filterChain.get<0>().coefficients->getMagnitudeForFrequency(pointToFreq, m_sampleRate);
        mag *= m_filterChain.get<1>().coefficients->getMagnitudeForFrequency(pointToFreq, m_sampleRate);
        mag *= m_filterChain.get<2>().coefficients->getMagnitudeForFrequency(pointToFreq, m_sampleRate);
        mag *= m_filterChain.get<3>().coefficients->getMagnitudeForFrequency(pointToFreq, m_sampleRate);

        difference += std::pow(m_drawnPoints[i] - mag, 2); // add squared difference
    }

    return difference;
}

void FilterSolver::updateCoefficientsBulk(const dlib::matrix<double, 3*NUM_FILTERS, 1>& variables) {
    const FilterParams fp0{ variables(0) , variables(1) , variables(2) };
    auto peak0Coefficients = makePeakFilter(fp0, m_sampleRate);
    updateCoefficients(m_filterChain.get<0>().coefficients, peak0Coefficients);

    const FilterParams fp1{ variables(3) , variables(4) , variables(5) };
    auto peak1Coefficients = makePeakFilter(fp1, m_sampleRate);
    updateCoefficients(m_filterChain.get<1>().coefficients, peak1Coefficients);

    const FilterParams fp2{ variables(6) , variables(7) , variables(8) };
    auto peak2Coefficients = makePeakFilter(fp2, m_sampleRate);
    updateCoefficients(m_filterChain.get<2>().coefficients, peak2Coefficients);

    const FilterParams fp3{ variables(9) , variables(10) , variables(11) };
    auto peak3Coefficients = makePeakFilter(fp3, m_sampleRate);
    updateCoefficients(m_filterChain.get<3>().coefficients, peak3Coefficients);
}

void FilterSolver::runSolver()
{
    // Initialize matrix
    dlib::matrix<double, 3*NUM_FILTERS, 1> variables;
    for (int fpos = 0; fpos < NUM_FILTERS; ++fpos) {
        float startingFreq = juce::mapToLog10(float(fpos + 1) / (NUM_FILTERS + 1), 20.f, 20000.f);
        variables(3 * fpos) = startingFreq; // Frequency guess
        variables(3 * fpos + 1) = 10.f; // Q guess
        variables(3 * fpos + 2) = 0.f; // Boost/Cut DB guess
    }

    updateCoefficientsBulk(variables);

    auto objectiveFunction = [this](const dlib::matrix<double, 3 * NUM_FILTERS, 1>& vars) {
        return this->solve(vars);
    };

    dlib::find_min_using_approximate_derivatives(
        dlib::bfgs_search_strategy(),
        dlib::objective_delta_stop_strategy(1e-3), // Play with this
        objectiveFunction, // Callable
        variables,
        -1, // Minimize
        5
    );

    for (int fpos = 0; fpos < NUM_FILTERS; ++fpos) {
        const FilterParams fp{ variables(3* fpos) , variables(3* fpos +1) , variables(3* fpos +2) };
        std::ostringstream oss;
        oss << fp;
        DBG("FilterParams: " << oss.str());
    }
}
