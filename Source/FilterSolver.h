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
#include <vector>
#include <Eigen/Core>

namespace fsolve {

constexpr int NUM_PARAMS = 3 * NUM_FILTERS;

struct FPSteps {
    float freqStep = 0.001f;
    float qStep = 0.001f;
    float boostCutDBStep = 0.001f;
};

using Eigen::VectorXf;

class FilterSolver {
public:
    FilterSolver(std::vector<float> drawnPoints, float sampleRate);
    ~FilterSolver();

    std::vector<FilterParams> runSolver();
    float operator()(const VectorXf& variables, VectorXf& grad);

    int iterationCount{ 0 };

private:
    void updateCoefficientsBulk(const VectorXf& variables);
    float getGrad(const VectorXf& vars, int pos, float stepSize);
    float solve(const VectorXf& variables);

    std::vector<float> m_drawnDecibels;
    float m_sampleRate;
    juce::dsp::ProcessorChain<Filter, Filter, Filter, Filter> m_filterChain; // TODO: Hardcoded to 4 for now
    FPSteps m_fpSteps;
};

}

