/*
  ==============================================================================

    FilterSolver.h
    Created: 30 Dec 2024 7:47:50pm
    Author:  Barak

  ==============================================================================
*/

#pragma once

#include "PluginProcessor.h"

#include <JuceHeader.h>
#include <dlib/optimization.h>
#include <vector>

class FilterSolver {
public:
    FilterSolver(std::vector<int> drawnPoints, double sampleRate);
    ~FilterSolver();

    double solve(const dlib::matrix<double, 3 * NUM_FILTERS, 1>& variables);

    void runSolver();

private:
    void updateCoefficientsBulk(const dlib::matrix<double, 3 * NUM_FILTERS, 1>& variables);

    std::vector<int> m_drawnPoints;
    double m_sampleRate;
    juce::dsp::ProcessorChain<Filter, Filter, Filter, Filter> m_filterChain; // TODO: Hardcoded to 4 for now
};