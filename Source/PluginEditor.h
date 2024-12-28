/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "PluginProcessor.h"

#include <vector>

class FilterCircle : public juce::Component {
public:
    FilterCircle();
    ~FilterCircle();
    void paint(juce::Graphics& g) override;
};

class EQGraphicComponent : public juce::Component {
public:
    EQGraphicComponent(PicassoEQAudioProcessor& ap);
    ~EQGraphicComponent();
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

    void updateFilterParamsFromCoords(int filterIndex, float eq_x, float eq_y);
    void updateResponseCurve();

private:
    juce::Path m_responseCurve;
    PicassoEQAudioProcessor& m_audioProcessor;
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
