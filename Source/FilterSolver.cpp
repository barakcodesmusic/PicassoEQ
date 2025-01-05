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

namespace fsolve {

using namespace LBFGSpp;

FilterSolver::FilterSolver(std::vector<int> drawnPoints, float sampleRate) :
    m_drawnPoints(drawnPoints),
    m_sampleRate(sampleRate)
{
}
FilterSolver::~FilterSolver() {}

float FilterSolver::solve(const VectorXf& variables) {

    updateCoefficientsBulk(variables);

    // Go through each drawn point, convert it to a frequency, get the filter response
    // and compare to expected response
    float difference = 0.f;
    for (int i = 0; i < m_drawnPoints.size(); ++i) {
        float mag = 1.f;
        float point = i / m_drawnPoints.size();
        float pointToFreq = juce::mapToLog10(point, FREQ_RANGE.first, FREQ_RANGE.second);

        mag *= m_filterChain.get<0>().coefficients->getMagnitudeForFrequency(pointToFreq, m_sampleRate);
        mag *= m_filterChain.get<1>().coefficients->getMagnitudeForFrequency(pointToFreq, m_sampleRate);
        mag *= m_filterChain.get<2>().coefficients->getMagnitudeForFrequency(pointToFreq, m_sampleRate);
        mag *= m_filterChain.get<3>().coefficients->getMagnitudeForFrequency(pointToFreq, m_sampleRate);

        difference += std::pow(m_drawnPoints[i] - mag, 2); // add squared difference
    }

    return difference;
}

float FilterSolver::getGrad(const VectorXf& vars, int pos, float stepSize) {
    VectorXf perturb_forw = vars;
    perturb_forw(pos) += stepSize;

    VectorXf perturb_back = vars;
    perturb_back(pos) -= stepSize;

    float solveForward = solve(perturb_forw);
    float solveBackward = solve(perturb_back);
    float grad = (solveForward - solveBackward) / (2 * stepSize);
    return grad;
}

float FilterSolver::operator()(const VectorXf& variables, VectorXf& grad) {

    for (int fpos = 0; fpos < NUM_PARAMS; fpos += 3) {
        // freq step
        grad[fpos] += getGrad(variables, fpos, this->m_fpSteps.freqStep);
        // q step
        grad[fpos + 1] = getGrad(variables, fpos + 1, this->m_fpSteps.qStep);
        // gain step
        grad[fpos + 2] = getGrad(variables, fpos + 2, this->m_fpSteps.boostCutDBStep);
    }
    // TODO: Update grad steps (for Q especially)
    //this->m_fpSteps.updateSteps(vars(i), vars(i + 1), vars(i + 2));

    return solve(variables);
}

void FilterSolver::updateCoefficientsBulk(const VectorXf& variables) {
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
    // Create starting vector
    VectorXf variables(NUM_PARAMS);
    variables.setZero();
    for (int fpos = 0; fpos < NUM_PARAMS; fpos+=3) {
        int filter = fpos / 3;
        float startingFreq = juce::mapToLog10(float(filter + 1) / (NUM_FILTERS + 1), 20.f, 20000.f);
        variables(fpos) = startingFreq; // Frequency guess
        variables(fpos + 1) = 1.f; // Q guess
        variables(fpos + 2) = 0.f; // Boost/Cut DB guess
    }

    updateCoefficientsBulk(variables);

    // Set up parameters
    LBFGSBParam<float> param;
    param.epsilon = 1e-6;
    param.max_iterations = 10;

    // Create solver and function object
    LBFGSBSolver<float> solver(param);

    VectorXf lb = VectorXf::Constant(NUM_PARAMS, 0.0f);
    VectorXf ub = VectorXf::Constant(NUM_PARAMS, 0.0f);
    // Set bounds for fcs
    for (int i = 0; i < NUM_FILTERS; i++)
    {
        int fpos = 3 * i;
        lb[fpos] = FREQ_RANGE.first;
        lb[fpos+1] = Q_RANGE.first;
        lb[fpos+2] = GAIN_RANGE.first;
        ub[fpos] = FREQ_RANGE.second;
        ub[fpos + 1] = Q_RANGE.second;
        ub[fpos + 2] = GAIN_RANGE.second;
    }

    // x will be overwritten to be the best point found
    float fx;
    int niter = solver.minimize(*this, variables, fx, lb, ub);

    for (int fpos = 0; fpos < NUM_FILTERS; ++fpos) {
        const FilterParams fp{ variables(3 * fpos) , variables(3 * fpos + 1) , variables(3 * fpos + 2) };
        std::ostringstream oss;
        oss << fp;
        DBG("FilterParams: " << oss.str());
    }
}

}