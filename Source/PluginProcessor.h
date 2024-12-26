/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>

#include "dsp/IRRFilter.h"
#include <unordered_map>


constexpr int NUM_FILTERS = 4;

template <int NumFilters, typename FilterType>
class FilterChain {
public:
    void reset(dsp::FilterParams& fp, int i, float sampleRate) {
        filters[i].reset(fp, sampleRate);
    }
    float processSample(float xn) {
        float mag = 1.f;
        for (auto& filter : filters) {
            mag *= filter.processSample(xn);
        }
        return mag;
    };

    std::array<FilterType, NumFilters> filters;
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

    dsp::FilterParams getUserFilterParams(int filterIndex);
    void updateFilters();

private:
    FilterChain<NUM_FILTERS, dsp::IRRFilter> leftChain;
    FilterChain<NUM_FILTERS, dsp::IRRFilter> rightChain;
    juce::dsp::Oscillator<float> osc;
    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (PicassoEQAudioProcessor)
};
