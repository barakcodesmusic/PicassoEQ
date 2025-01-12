/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "PluginProcessor.h"

#include <vector>

using EQPoint = juce::Point<int>;

template <int NumFilters, typename Component, typename FilterChain>
struct FilterInterface {
    // Returns >0 if x/y inside bounds of one of the circles
    // this FilterInfo object owns
    int inBoundsOfCircle(const EQPoint& p) {
        for (int i = 0; i < NumFilters; ++i) {
            if (fcomps[i].getBounds().contains(p)) {
                return i;
            }
        }
        return -1;
    };
    void setCoeffs(FilterParams& fp, int i, float sampleRate) {
        fchain.setCoeffs(fp, i, sampleRate);
    }
    int getFilterIndex(Component& comp) {
        for (int i = 0; i < NumFilters; ++i) {
            if (fcomps[i] == comp) {
                return i;
            }
        }
        jassert(false); // Should never reach here, means we didn't find passed down filter
    }

    int size = NumFilters;
    FilterChain fchain;
    std::array<Component, NumFilters> fcomps;
};

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
    void setFilterAnchorPositions();
    void resized() override;
    void mouseDrag(const juce::MouseEvent& event) override;
    void mouseDown(const juce::MouseEvent& event) override;
    void mouseUp(const juce::MouseEvent& event) override;
    void mouseDoubleClick(const juce::MouseEvent& event) override;

    void drawTextLabels(juce::Graphics& g);
    
    void drawBackgroundGrid(juce::Graphics& g);

    std::vector<float> getXs(const std::vector<float>& freqs, float left, float width);
    juce::Rectangle<int> getRenderArea();
    juce::Rectangle<int> getAnalysisArea();

    std::vector<float> getFrequencies();
    std::vector<float> getGains();

    void updateChain();

private:
    void updateDrawnCurve();
    void updateResponseCurve();
    void adjustFiltersAtClickPoint(const EQPoint& filterPos);
    void updateFilterParamsFromCoords(int filterIndex, const EQPoint& filterPos);
    void resetCurveDraw(const EQPoint& mouseDownPoint);

    PicassoEQAudioProcessor& m_audioProcessor;
    juce::Path m_responseCurve;
    juce::Path m_drawCurve;
    bool m_drawing;
    std::vector<int> m_drawnPoints;
    int m_prevX = -1; // TODO: Clean up, should be in seperate class probably
    bool m_drewToAxis = false; // TODO: Clean up, should be in seperate class probably
    FilterInterface<NUM_FILTERS, FilterCircle, MonoChain> m_filterInterface;
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
