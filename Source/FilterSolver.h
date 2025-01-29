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
#include <Eigen/Core>
#include <vector>
#include <functional>

namespace fsolve {

constexpr int PARAMS_PER_FILTER = 3;
constexpr int NUM_PARAMS = PARAMS_PER_FILTER * NUM_FILTERS;

struct FPSteps {
    float freqStep = 0.00001f;
    float qStep = 0.0001f;
    float boostCutDBStep = 0.0001f;
};

using Eigen::VectorXf;

using CoefficientSolver = std::function<Coefficients(const FilterParams&, double)>;
using Callback = std::function<void(const int, const FilterParams)>;

class FilterSolverThread : public juce::Thread {
public:
    FilterSolverThread(
        std::vector<float>&& dbsToSolve,
        CoefficientSolver&& coefficientSolver, 
        FilterParams&& startingParams,
        const int filterIndex,
        const double sampleRate,
        Callback cb);
    virtual ~FilterSolverThread();

    virtual void run();

    float operator()(const VectorXf& variables, VectorXf& grad);

    int iterationCount{ 0 };

    FilterParams getFilterParams() const { return m_filterParams; };

private:
    void solverUpdateCoefficients(const VectorXf& variables);
    float getGrad(const VectorXf& vars, int pos, float stepSize);
    float solve(const VectorXf& variables);

    std::vector<float> m_dbsToSolve;
    Filter m_filter;
    CoefficientSolver m_coefficientSolver;
    FilterParams m_filterParams;
    const int m_filterIndex;
    const double m_sampleRate;
    Callback m_cb;
    FPSteps m_fpSteps;
};

}

