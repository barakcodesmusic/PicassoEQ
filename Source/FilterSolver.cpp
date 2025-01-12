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

FilterSolver::FilterSolver(std::vector<float> drawnPoints, float sampleRate) :
    m_drawnDecibels(drawnPoints),
    m_sampleRate(sampleRate)
{
}
FilterSolver::~FilterSolver() {}

float FilterSolver::solve(const VectorXf& variables) {

    updateCoefficientsBulk(variables);

    // Go through each drawn point, convert it to a frequency, get the filter response
    // and compare to expected response
    float loss = 0.f;
    for (int i = 0; i < m_drawnDecibels.size(); ++i) {
        float mag = 1.f;
        float point = float(i) / m_drawnDecibels.size();
        float pointToFreq = juce::mapToLog10(point, FREQ_RANGE.first, FREQ_RANGE.second);

        mag *= m_filterChain.get<0>().coefficients->getMagnitudeForFrequency(pointToFreq, m_sampleRate);
        mag *= m_filterChain.get<1>().coefficients->getMagnitudeForFrequency(pointToFreq, m_sampleRate);
        mag *= m_filterChain.get<2>().coefficients->getMagnitudeForFrequency(pointToFreq, m_sampleRate);
        mag *= m_filterChain.get<3>().coefficients->getMagnitudeForFrequency(pointToFreq, m_sampleRate);
        
        float calcDecibels = juce::Decibels::gainToDecibels(mag);
        loss += std::abs(m_drawnDecibels[i] - calcDecibels); // absolute difference
    }

    return loss;
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

    if (iterationCount % 5 == 0) {
        for (int fpos = 0; fpos < NUM_PARAMS; fpos += 3) {
            // freq step
            grad[fpos] = getGrad(variables, fpos, this->m_fpSteps.freqStep);
            // q step
            grad[fpos + 1] = getGrad(variables, fpos + 1, this->m_fpSteps.qStep);
            // gain step
            grad[fpos + 2] = getGrad(variables, fpos + 2, this->m_fpSteps.boostCutDBStep);
        }
    }

    // TODO: Update grad steps (for Q especially)
    //this->m_fpSteps.updateSteps(vars(i), vars(i + 1), vars(i + 2));
    iterationCount += 1;
    return solve(variables);
}

void FilterSolver::updateCoefficientsBulk(const VectorXf& variables) {
    // TODO: Only make when necessary
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

std::vector<FilterParams> FilterSolver::runSolver()
{
    // Create starting vector
    VectorXf variables(NUM_PARAMS);
    variables.setZero();
    for (int fpos = 0; fpos < NUM_PARAMS; fpos+=3) {
        int filter = fpos / 3;
        float startingFreq = juce::mapToLog10(float(filter + 1) / (NUM_FILTERS + 1), 20.f, 20000.f);
        variables(fpos) = startingFreq; // Frequency guess
        variables(fpos + 1) = 2.f; // Q guess
        variables(fpos + 2) = 0.f; // Boost/Cut DB guess
    }

    updateCoefficientsBulk(variables);

    // Set up parameters
    LBFGSBParam<float> param;
    param.epsilon = 1e-4;
    param.epsilon_rel = 1e-4;
    param.max_iterations = 150;
    param.max_linesearch = 1000;
    param.ftol = 1e-2;
    param.wolfe = 0.7;
    //param.ftol = 0.4;
    //param.wolfe = 0.9;

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
    try {
        float fx;
        int niter = solver.minimize(*this, variables, fx, lb, ub);
    }
    catch (...) {
        DBG("Solver quit");
    }
    

    std::vector<FilterParams> out;
    for (int fpos = 0; fpos < NUM_FILTERS; ++fpos) {
        FilterParams fp{ variables(3 * fpos) , variables(3 * fpos + 1) , variables(3 * fpos + 2) };
        out.push_back(fp);
        std::ostringstream oss;
        oss << out[fpos];
        DBG("FilterParams: " << oss.str());
    }
    return out;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

FilterSolverThread::FilterSolverThread(
    std::vector<float>&& dbsToSolve,
    CoefficientSolver&& coefficientSolver, 
    FilterParams&& startingParams,
    PicassoEQAudioProcessor& audioProcessor,
    const juce::String& filterName) :
    juce::Thread(filterName),
    m_dbsToSolve(std::move(dbsToSolve)),
    m_coefficientSolver(std::move(coefficientSolver)),
    m_filterParams(std::move(startingParams)),
    m_audioProcessor(audioProcessor),
    m_filterName(filterName)
{
}

FilterSolverThread::~FilterSolverThread()
{
}

void FilterSolverThread::solverUpdateCoefficients(const VectorXf& variables) {
    // TODO: Only make when necessary
    const FilterParams guessFP{ variables(0) , variables(1) , variables(2) };
    auto peakCoefficients = m_coefficientSolver(guessFP, m_audioProcessor.getSampleRate());
    updateCoefficients(m_filter.coefficients, peakCoefficients);
}

float FilterSolverThread::getGrad(const VectorXf& vars, int pos, float stepSize) {
    VectorXf perturb_forw = vars;
    perturb_forw(pos) += stepSize;

    VectorXf perturb_back = vars;
    perturb_back(pos) -= stepSize;

    float solveForward = solve(perturb_forw);
    float solveBackward = solve(perturb_back);
    float grad = (solveForward - solveBackward) / (2 * stepSize);
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
        float mag = m_filter.coefficients->getMagnitudeForFrequency(pointToFreq, m_audioProcessor.getSampleRate());
        float calcDecibels = juce::Decibels::gainToDecibels(mag);
        loss += std::abs(m_dbsToSolve[i] - calcDecibels); // absolute difference
    }
    return loss;
}

float FilterSolverThread::operator()(const VectorXf& variables, VectorXf& grad) {

    if (iterationCount % 5 == 0) {
        grad[0] = getGrad(variables, 0, this->m_fpSteps.freqStep);
        // q step
        grad[1] = getGrad(variables, 1, this->m_fpSteps.qStep);
        // gain step
        grad[2] = getGrad(variables, 2, this->m_fpSteps.boostCutDBStep);   
    }

    // TODO: Update grad steps (for Q especially)
    //this->m_fpSteps.updateSteps(vars(i), vars(i + 1), vars(i + 2));
    iterationCount += 1;
    return solve(variables);
}

void FilterSolverThread::run()
{
    // Create starting vector
    VectorXf variables(PARAMS_PER_FILTER);
    variables.setZero();

    variables(0) = m_filterParams.cutoffFreq; // Frequency guess
    variables(1) = m_filterParams.q; // Q guess
    variables(2) = m_filterParams.boostCutDB; // Boost/Cut DB guess
   
    solverUpdateCoefficients(variables);

    // Set up parameters
    LBFGSBParam<float> param;
    param.epsilon = 1e-4;
    param.epsilon_rel = 1e-4;
    param.max_iterations = 250;
    param.max_linesearch = 1000;
    param.ftol = 1e-2;
    param.wolfe = 0.7;
    //param.ftol = 0.4;
    //param.wolfe = 0.9;

    // Create solver and function object
    LBFGSBSolver<float> solver(param);

    VectorXf lb = VectorXf::Constant(PARAMS_PER_FILTER, 0.0f);
    VectorXf ub = VectorXf::Constant(PARAMS_PER_FILTER, 0.0f);

    // Set bounds for fcs
    lb[0] = FREQ_RANGE.first;
    lb[1] = Q_RANGE.first;
    lb[2] = GAIN_RANGE.first;
    ub[0] = FREQ_RANGE.second;
    ub[1] = Q_RANGE.second;
    ub[2] = GAIN_RANGE.second;
    
    // x will be overwritten to be the best point found
    try {
        float fx;
        int niter = solver.minimize(*this, variables, fx, lb, ub);
    }
    catch (...) {
        DBG("Solver quit");
    }

    m_filterParams.cutoffFreq = variables(0);
    m_filterParams.q = variables(1);
    m_filterParams.boostCutDB = variables(2);

    auto mapParamToFrac = [](auto param, auto start, auto end) {
        return juce::jmap(param, start, end, 0.f, 1.f);
    };

    auto mapFreqToFrac = [&](auto freq) {return mapParamToFrac(freq, FREQ_RANGE.first, FREQ_RANGE.second);};
    auto mapQToFrac = [&](auto q) {return mapParamToFrac(q, Q_RANGE.first, Q_RANGE.second);};
    auto mapGainToFrac = [&](auto gain) {return mapParamToFrac(gain, GAIN_RANGE.first, GAIN_RANGE.second);};

    m_audioProcessor.apvts.getParameter(m_filterName + "LowCut Freq")->setValueNotifyingHost(mapFreqToFrac(m_filterParams.cutoffFreq));
    m_audioProcessor.apvts.getParameter(m_filterName + "Q")->setValueNotifyingHost(mapQToFrac(m_filterParams.q));
    m_audioProcessor.apvts.getParameter(m_filterName + "BoostCutDB")->setValueNotifyingHost(mapGainToFrac(m_filterParams.boostCutDB));

    std::ostringstream oss;
    oss << m_filterParams;
    DBG("FilterParams: " << oss.str());
}

}