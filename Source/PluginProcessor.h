/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>

#include "TypeHelper.h"

constexpr int NUM_FILTERS = 4;
constexpr std::pair<float, float> FREQ_RANGE{ 20.f, 20000.f };
constexpr std::pair<float, float> Q_RANGE{ 0.1f, 18.f };
constexpr std::pair<float, float> GAIN_RANGE{ -24.f, 24.f };

auto mapFracToParam = [](auto frac, auto start, auto end) {
    return juce::jmap(frac, start, end);
};
auto mapParamToFrac = [](auto param, auto start, auto end) {
    return juce::jmap(param, start, end, 0.f, 1.f);
};

auto mapFreqToFrac = [&](auto freq) {return mapParamToFrac(freq, FREQ_RANGE.first, FREQ_RANGE.second);};
auto mapQToFrac = [&](auto q) {return mapParamToFrac(q, Q_RANGE.first, Q_RANGE.second);};
auto mapDBToFrac = [&](auto gain) {return mapParamToFrac(gain, GAIN_RANGE.first, GAIN_RANGE.second);};

auto mapFracToFreq = [&](auto frac) {return mapFracToParam(frac, FREQ_RANGE.first, FREQ_RANGE.second);};
auto mapFracToQ = [&](auto frac) {return mapFracToParam(frac, Q_RANGE.first, Q_RANGE.second);};
auto mapFracToDB = [&](auto frac) {return mapFracToParam(frac, GAIN_RANGE.first, GAIN_RANGE.second);};

const juce::String CUTOFF("Cutoff");
const juce::String BOOSTCUTDB("BoostCutDB");
const juce::String Q("Q");

extern juce::String cutoffParamFromIndex(int i);
extern juce::String gainDBParamFromIndex(int i);
extern juce::String qParamFromIndex(int i);

enum Slope
{
    Slope_12,
    Slope_24,
    Slope_36,
    Slope_48
};

const Slope SLOPE = Slope_24; // TODO: Hardcoded order!

using Filter = juce::dsp::IIR::Filter<float>;
using CutFilter = juce::dsp::ProcessorChain<Filter, Filter, Filter, Filter>;
using MonoChain = juce::dsp::ProcessorChain<Filter, Filter, Filter, Filter>;

struct FilterParams {
    float cutoffFreq;
    float q;
    float boostCutDB;
    //FilterAlgorithm fa;

    bool operator==(const FilterParams& r) const {
        return std::tie(cutoffFreq, q, boostCutDB) == std::tie(r.cutoffFreq, r.q, r.boostCutDB);
        //return std::tie(cutoffFreq, q, boostCutDB, fa) == std::tie(r.cutoffFreq, r.q, r.boostCutDB, r.fa);
    }
    bool operator!=(const FilterParams& r) const {
        return !(*this == r);
    }

    friend std::ostream& operator<<(std::ostream& os, const FilterParams& obj) {
        os << "Cutoff: " << obj.cutoffFreq << ", Q: " << obj.q << " Boost/Cut DB: " << obj.boostCutDB;
        return os;
    }
};

using Coefficients = Filter::CoefficientsPtr;
void updateCoefficients(Coefficients& old, const Coefficients& replacements);

Coefficients makePeakFilter(const FilterParams& filterParams, double sampleRate);

template<int Index, typename ChainType, typename CoefficientType>
void update(ChainType& chain, const CoefficientType& coefficients)
{
    updateCoefficients(chain.template get<Index>().coefficients, coefficients[Index]);
    chain.template setBypassed<Index>(false);
}

template<typename ChainType, typename CoefficientType>
void updateCutFilter(ChainType& chain,
    const CoefficientType& coefficients,
    const Slope& slope)
{
    chain.template setBypassed<0>(true);
    chain.template setBypassed<1>(true);
    chain.template setBypassed<2>(true);
    chain.template setBypassed<3>(true);

    switch (slope)
    {
    case Slope_48:
    {
        update<3>(chain, coefficients);
    }
    case Slope_36:
    {
        update<2>(chain, coefficients);
    }
    case Slope_24:
    {
        update<1>(chain, coefficients);
    }
    case Slope_12:
    {
        update<0>(chain, coefficients);
    }
    }
}

inline auto makeLowCutFilter(const FilterParams& filterParams, double sampleRate)
{
    return juce::dsp::FilterDesign<float>::designIIRHighpassHighOrderButterworthMethod(filterParams.cutoffFreq,
        sampleRate,
        2 * (SLOPE + 1)); // TODO: Hardcoded order!
}

inline auto makeHighCutFilter(const FilterParams& filterParams, double sampleRate)
{
    return juce::dsp::FilterDesign<float>::designIIRLowpassHighOrderButterworthMethod(filterParams.cutoffFreq,
        sampleRate,
        2 * (SLOPE + 1)); // TODO: Hardcoded order!
}

enum ChainPositions
{
    LowCut,
    Peak1,
    Peak2,
    HighCut
};

class PicassoEQAudioProcessor  : public juce::AudioProcessor
{
public:
    //==============================================================================
    PicassoEQAudioProcessor();
    ~PicassoEQAudioProcessor() override;

    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

   #ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
   #endif

    void processBlock (juce::AudioBuffer<float>&, juce::MidiBuffer&) override;

    //==============================================================================
    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const juce::String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const juce::String getProgramName (int index) override;
    void changeProgramName (int index, const juce::String& newName) override;

    //==============================================================================
    void getStateInformation (juce::MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;

    static juce::AudioProcessorValueTreeState::ParameterLayout createParameterLayout();
    juce::AudioProcessorValueTreeState apvts{ *this, nullptr, "Parameters", createParameterLayout() };

    FilterParams getUserFilterParams(int filterIndex) const;

    template <ChainPositions position>
    void updatePeakFilter(const FilterParams& filterParams)
    {
        auto peakCoefficients = makePeakFilter(filterParams, getSampleRate());

        //leftChain.setBypassed<ChainPositions::Peak>(filterParams.peakBypassed);
        //rightChain.setBypassed<ChainPositions::Peak>(filterParams.peakBypassed);

        updateCoefficients(leftChain.get<position>().coefficients, peakCoefficients);
        updateCoefficients(rightChain.get<position>().coefficients, peakCoefficients);
    }
    void updateLowCutFilters(const FilterParams& filterParams);
    void updateHighCutFilters(const FilterParams& filterParams);
    void updateFilters();

private:
    MonoChain leftChain, rightChain;
    juce::dsp::Oscillator<float> osc;
    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (PicassoEQAudioProcessor)
};
