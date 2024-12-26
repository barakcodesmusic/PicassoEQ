/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "PluginProcessor.h"

#include <vector>

using APVTS = juce::AudioProcessorValueTreeState;

struct FilterCircle : juce::Component {
    FilterCircle();
    ~FilterCircle();
    void paint(juce::Graphics& g) override;
};

class EQGraphicComponent : public juce::Component {
public:
    EQGraphicComponent(PicassoEQAudioProcessor& ap);
    ~EQGraphicComponent();
    std::vector<juce::Component*> getComps();
    void paint(juce::Graphics& g) override;
    void resized() override;
    void mouseDrag(const juce::MouseEvent& event) override;

    void drawTextLabels(juce::Graphics& g);
    void drawBackgroundGrid(juce::Graphics& g);

    std::vector<float> getXs(const std::vector<float>& freqs, float left, float width);
    juce::Rectangle<int> getRenderArea();
    juce::Rectangle<int> getAnalysisArea();

    std::vector<float> getFrequencies();
    std::vector<float> getGains();

    void setNewFilterParams(float eq_x, float eq_y);
    void updateResponseCurve();

private:
    juce::Path m_responseCurve;
    std::vector<FilterCircle> m_filterCircles;
    PicassoEQAudioProcessor& m_audioProcessor;

    dsp::IRRFilter m_lowPassFilter;

};

//==============================================================================
/**
*/
class PicassoEQAudioProcessorEditor  : public juce::AudioProcessorEditor
{
public:
    PicassoEQAudioProcessorEditor (PicassoEQAudioProcessor&);
    ~PicassoEQAudioProcessorEditor() override;

    //==============================================================================
    void paint (juce::Graphics&) override;
    void resized() override;

private:
    // This reference is provided as a quick way for your editor to
    // access the processor object that created it.
    PicassoEQAudioProcessor& audioProcessor;

    EQGraphicComponent eqGraphicComponent;

    std::vector<juce::Component*> getComps();

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (PicassoEQAudioProcessorEditor)
};
