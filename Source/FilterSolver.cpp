/*
  ==============================================================================

    FilterSolver.cpp
    Created: 30 Dec 2024 7:47:50pm
    Author:  Barak

  ==============================================================================
*/

#include "FilterSolver.h"

#include <iostream>
#include <LBFGSB.h>

#include <stdexcept>

namespace fsolve {

using namespace LBFGSpp;

FilterSolverThread::FilterSolverThread(
    std::vector<float>&& dbsToSolve,
    CoefficientSolver&& coefficientSolver, 
    FilterParams&& startingParams,
    const int filterIndex,
    const double sampleRate,
    Callback callback) :
    juce::Thread(std::to_string(filterIndex)),
    m_dbsToSolve(std::move(dbsToSolve)),
    m_coefficientSolver(std::move(coefficientSolver)),
    m_filterParams(std::move(startingParams)),
    m_filterIndex(filterIndex),
    m_sampleRate(sampleRate),
    m_cb(callback)
{
}

FilterSolverThread::~FilterSolverThread()
{
}

void FilterSolverThread::solverUpdateCoefficients(const VectorXf& variables) {
    // TODO: Only make when necessary
    const FilterParams guessFP{ mapFracToFreq(variables(0)) , mapFracToQ(variables(1)) , mapFracToDB(variables(2)) };
    auto peakCoefficients = m_coefficientSolver(guessFP, m_sampleRate);
    updateCoefficients(m_filter.coefficients, peakCoefficients);
}

float FilterSolverThread::getGrad(const VectorXf& vars, int pos, float stepSize) {
    VectorXf perturb_forw = vars;
    perturb_forw(pos) += stepSize;

    VectorXf perturb_back = vars;
    perturb_back(pos) -= stepSize;

    float f = solve(perturb_forw);
    float b = solve(perturb_back);

    float grad = (f - b) / (2.f * stepSize);
    return grad;
}

float FilterSolverThread::solve(const VectorXf& variables) {

    solverUpdateCoefficients(variables);

    // Go through each drawn point, convert it to a frequency, get the filter response
    // and compare to expected response
    float loss = 0.f;
    for (int i = 0; i < m_dbsToSolve.size(); ++i) {
        float point = float(i) / m_dbsToSolve.size();
        float pointToFreq = juce::mapToLog10(point, FREQ_RANGE.first, FREQ_RANGE.second);
        float mag = m_filter.coefficients->getMagnitudeForFrequency(pointToFreq, m_sampleRate);
        float calcDecibels = juce::Decibels::gainToDecibels(mag);
        float dist = mapDBToFrac(m_dbsToSolve[i]) - mapDBToFrac(calcDecibels);
        loss += std::abs(dist); // absolute difference
    }
    return loss;
}

float FilterSolverThread::operator()(const VectorXf& variables, VectorXf& grad) {

    // freq step
    grad[0] = getGrad(variables, 0, this->m_fpSteps.freqStep);
    // q step
    grad[1] = getGrad(variables, 1, this->m_fpSteps.qStep);
    // gain step
    grad[2] = getGrad(variables, 2, this->m_fpSteps.boostCutDBStep);   

    iterationCount += 1;
    return solve(variables);
}

void FilterSolverThread::run()
{
    // Create starting vector
    VectorXf variables(PARAMS_PER_FILTER);
    variables.setZero();

    variables(0) = mapFreqToFrac(m_filterParams.cutoffFreq); // Frequency guess
    variables(1) = mapQToFrac(m_filterParams.q); // Q guess
    variables(2) = mapDBToFrac(m_filterParams.boostCutDB); // Boost/Cut DB guess
   
    solverUpdateCoefficients(variables);

    // Set up parameters
    LBFGSBParam<float> param;
    param.epsilon = 1e-10;
    //param.delta = 1e-10;
    //param.past = 20;
    //param.epsilon_rel = 1e-4;
    param.max_iterations = 200;
    param.max_linesearch = 500;
    param.m = 20;
    //param.ftol = 0.1;
    param.min_step = 1e-5;
    //param.ftol = 1e-2;
    //param.wolfe = 0.7;

    VectorXf lb = VectorXf::Constant(PARAMS_PER_FILTER, 0.0f);
    VectorXf ub = VectorXf::Constant(PARAMS_PER_FILTER, 0.0f);

    // Set bounds for fcs
    lb[0] = 0;
    lb[1] = 0;
    lb[2] = 0;
    ub[0] = 1;
    ub[1] = 1;
    ub[2] = 1;
    
    try {
        // Create solver and function object
        LBFGSBSolver<float> solver(param);
        float fx = solve(variables);
        int result = solver.minimize(*this, variables, fx, lb, ub);
    }
    catch (const std::exception& e) {
        DBG("Solver quit");
    }

    m_filterParams.cutoffFreq = mapFracToFreq(variables(0));
    m_filterParams.q = mapFracToQ(variables(1));
    m_filterParams.boostCutDB = mapFracToDB(variables(2));

    m_cb(m_filterIndex, m_filterParams);

    std::ostringstream oss;
    oss << m_filterParams;
    DBG("FilterParams: " << oss.str());
}

}