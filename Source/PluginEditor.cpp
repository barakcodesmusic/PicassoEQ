/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
PicassoEQAudioProcessorEditor::PicassoEQAudioProcessorEditor (PicassoEQAudioProcessor& p)
    : AudioProcessorEditor (&p), audioProcessor (p), 
    eqGraphicComponent(audioProcessor)
{
    // Make sure that before the constructor has finished, you've set the
    // editor's size to whatever you need it to be.
    setSize (1280, 720);

    for (auto comp : getComps()) {
        addAndMakeVisible(comp);
    }
}

PicassoEQAudioProcessorEditor::~PicassoEQAudioProcessorEditor()
{
}

//==============================================================================
void PicassoEQAudioProcessorEditor::paint (juce::Graphics& g)
{
    using namespace juce;

    g.fillAll(Colours::blue);
    g.setFont (juce::FontOptions (15.0f));

    auto bounds = getLocalBounds();
    auto center = bounds.getCentre();

    String title = { "PicassoEQ" };
    g.setFont(25);
    auto titleWidth = g.getCurrentFont().getStringWidth(title);

    g.setColour(Colour(255u, 154u, 1u));
    g.drawFittedText(title, bounds, juce::Justification::centredTop, 1);
}

void PicassoEQAudioProcessorEditor::resized()
{
    auto bounds = getLocalBounds();

    auto eqGraphicArea = bounds;
    eqGraphicArea.removeFromTop(40); // JUCE_LIVE_CONSTANT
    eqGraphicArea.removeFromBottom(40);
    eqGraphicArea.removeFromBottom(30);
    eqGraphicArea.removeFromBottom(30);

    eqGraphicComponent.setBounds(eqGraphicArea);
}

std::vector<juce::Component*> PicassoEQAudioProcessorEditor::getComps()
{
    return {
        &eqGraphicComponent
    };
}

EQGraphicComponent::EQGraphicComponent(PicassoEQAudioProcessor& ap) :
    m_audioProcessor(ap)
{
    jassert(m_filterInterface.size >= 2);
    
    // TODO: Make this tied to NUM_FILTERS parameter?
    std::array<float, NUM_FILTERS> startFreqs{40, 300, 2000, 12000};
    // std::array<float, NUM_FILTERS> startFreqs{ 40 };
    for (int i = 0; i < m_filterInterface.size; ++i) {
        // TODO: check if default value exists?
        std::string filterID{ "Filter" + std::to_string(i) };
        if (i == 0) {
            float normalized = (dsp::FilterAlgorithm::kHPF2) / float(dsp::stringToFilterAlgorithm.size());
            m_audioProcessor.apvts.getParameter(filterID + "Filter Algorithm")->setValueNotifyingHost(normalized);
        }
        else if (i == m_filterInterface.size - 1) {
            float normalized = dsp::FilterAlgorithm::kLPF2 / float(dsp::stringToFilterAlgorithm.size());
            m_audioProcessor.apvts.getParameter(filterID + "Filter Algorithm")->setValueNotifyingHost(normalized);
        }
        else {
            float normalized = (dsp::FilterAlgorithm::kCQParaEQ + 1) / float(dsp::stringToFilterAlgorithm.size());
            m_audioProcessor.apvts.getParameter(filterID + "Filter Algorithm")->setValueNotifyingHost(normalized);
        }
        
        //m_audioProcessor.apvts.getParameter(filterID + "LowCut Freq")->setValueNotifyingHost(startFreqs[i] / 20000.f);
        //m_audioProcessor.apvts.getParameter(filterID + "Q")->setValueNotifyingHost(0.5f);   
        //m_audioProcessor.apvts.getParameter(filterID + "BoostCutDB")->setValueNotifyingHost(0.45f);
    }

    for (auto& comp : m_filterInterface.fcomps) {
        addAndMakeVisible(comp);
    }
    repaint();
}

EQGraphicComponent::~EQGraphicComponent()
{
}

void EQGraphicComponent::drawTextLabels(juce::Graphics& g)
{
    using namespace juce;
    g.setColour(Colours::lightgrey);
    const int fontHeight = 10;
    g.setFont(fontHeight);

    auto renderArea = getAnalysisArea();
    auto left = renderArea.getX();
    auto top = renderArea.getY();
    auto bottom = renderArea.getBottom();
    auto width = renderArea.getWidth();

    auto freqs = getFrequencies();
    auto xs = getXs(freqs, left, width);

    for (int i = 0; i < freqs.size(); ++i)
    {
        auto f = freqs[i];
        auto x = xs[i];

        bool addK = false;
        String str;
        if (f > 999.f)
        {
            addK = true;
            f /= 1000.f;
        }

        str << f;
        if (addK)
            str << "k";
        str << "Hz";

        auto textWidth = g.getCurrentFont().getStringWidth(str);

        Rectangle<int> r;

        r.setSize(textWidth, fontHeight);
        r.setCentre(x, 0);
        r.setY(1);

        g.drawFittedText(str, r, juce::Justification::centred, 1);
    }

    auto gain = getGains();

    for (auto gDb : gain)
    {
        auto y = jmap(gDb, -24.f, 24.f, float(bottom), float(top));

        // DB
        String str;
        if (gDb > 0)
            str << "+";
        str << gDb;

        auto textWidth = g.getCurrentFont().getStringWidth(str);

        Rectangle<int> r;
        r.setSize(textWidth, fontHeight);
        r.setX(1);
        r.setCentre(r.getCentreX(), y);

        g.setColour(gDb == 0.f ? Colour(0u, 172u, 1u) : Colours::lightgrey);

        g.drawFittedText(str, r, juce::Justification::centredLeft, 1);
    }
}

void EQGraphicComponent::paint(juce::Graphics& g)
{
    using namespace juce;

    auto bounds = getLocalBounds();

    g.fillAll(Colours::black);

    drawBackgroundGrid(g);

    drawTextLabels(g);

    g.setColour(Colours::white);
    g.strokePath(m_responseCurve, PathStrokeType(2.f));

    updateResponseCurve();
}

void EQGraphicComponent::resized() {
    auto bounds = getAnalysisArea();

    auto mapFPsToPositions = [&bounds](const auto& fp){

        //DBG(width << " " << top << " " << bottom << " " << centerY);

        float midPointY = bounds.getCentreY();
        float bottom = bounds.getBottom();
        float top = bounds.getY();
        
        float posX = bounds.getX() + bounds.getWidth() * juce::mapFromLog10(fp.cutoffFreq, 20.f, 20000.f);
        float posY = 0.0;
        if (dsp::IRRFilter::algorithmRequiresGain(fp.fa)) {
            // TODO: Figure out how to calculate q?
            posY = juce::jmap(fp.boostCutDB, -24.f, 24.f, bottom, top);
        }
        else {
            if (fp.q > 1) { // TODO: Not accurate midpoint?
                posY = juce::jmap(fp.q, 1.f, 18.f, midPointY, top); // 1 -> 18
            }
            else {
                posY = juce::jmap(fp.q, 0.1f, 1.f, bottom, midPointY); // .1 -> 1
            }
        }
        
        //DBG("Posx: " << posX << " Posy: " << posY << " width: " << bounds.getWidth() << " height: " << bounds.getHeight() << " top: " << bounds.getY() << " bottom " << bounds.getBottom());
        
        return juce::Rectangle<int>{static_cast<int>(posX), static_cast<int>(posY), 25, 25};
    };

    for (int i = 0; i < m_filterInterface.size; ++i) {
        auto& fp = m_audioProcessor.getUserFilterParams(i);
        //DBG("index: " << i << " cutoff: " << fp.cutoffFreq << " q: " << fp.q << " fa: " << fp.fa << " boost/cut: " << fp.boostCutDB);
        juce::Rectangle<int> filterBound = mapFPsToPositions(fp);
        m_filterInterface.fcomps[i].setBounds(filterBound);
    }
    repaint();
}

void EQGraphicComponent::updateFilterParamsFromCoords(int filterIndex, float eq_x, float eq_y) {
    // TODO: Make this more generic. For peak filters need to pass down gain
    auto bounds = getAnalysisArea();

    float height = bounds.getHeight();
    float width = bounds.getWidth();
    float top = bounds.getY();
    float bottom = bounds.getBottom();

    float cutoffFreq = juce::mapToLog10(eq_x/width, 20.f/20000.f, 1.f);
    float boostCutDB = juce::jmap(eq_y, bottom, top, 0.f, 1.f);

    auto& fp = m_audioProcessor.getUserFilterParams(filterIndex);

    float q = 0.0;
    float center = (top + bottom) / 2;
    if (eq_y < center) { // In upper half
        q = juce::jmap(eq_y, center, top, 1.f / 18.f, 1.f); // 1 -> 18
    }
    else {
        q = juce::jmap(eq_y, bottom, center, 1.f/180.f, 1.f / 18.f); // .1 -> 1
    }

    // Update Audio Processor // TODO: Maybe FilterInterface has a reference to AudioProcessor?
    std::string filterID{ "Filter" + std::to_string(filterIndex) };
    m_audioProcessor.apvts.getParameter(filterID + "LowCut Freq")->setValueNotifyingHost(cutoffFreq);

    if (m_audioProcessor.getUserFilterParams(filterIndex).fa == dsp::FilterAlgorithm::kCQParaEQ) {
        q = 1.f / 18.f;
    }
    m_audioProcessor.apvts.getParameter(filterID + "Q")->setValueNotifyingHost(q);
    m_audioProcessor.apvts.getParameter(filterID + "BoostCutDB")->setValueNotifyingHost(boostCutDB);

    //DBG("Cutoff: " << cutoffFreq << ", Q : " << q << ", Boost/Cut : " << boostCutDB);
}

void EQGraphicComponent::updateResponseCurve()
{
    using namespace juce;
    auto responseArea = getAnalysisArea();

    auto w = responseArea.getWidth();

    std::vector<double> mags;

    mags.resize(w);

    // Set UI side filter from audioProcessor.filterParams
    for (size_t i = 0; i < m_filterInterface.size; ++i) {
        auto& fp = m_audioProcessor.getUserFilterParams(i);
        m_filterInterface.setCoeffs(fp, i, m_audioProcessor.getSampleRate());
    }
    
    
    for (int i = 0; i < w; ++i) {
        float xPos = float(i) / float(w);
        auto freq = mapToLog10(xPos, 20.f, 20000.f);
        float mag = m_filterInterface.fchain.getMagnitudeForFrequency(freq, m_audioProcessor.getSampleRate());
        mags[i] = Decibels::gainToDecibels(mag);
        //DBG("Mag at freq: " << std::to_string(freq) << ": " << std::to_string(mags[i]));
    }

    const float outputMinY = responseArea.getBottom();
    const float outputMaxY = responseArea.getY();
    //DBG("Bottom: " << outputMinY << " Top: " << outputMaxY << " Center: " << responseArea.getCentreY() << " Height: " << responseArea.getHeight());
    auto map = [outputMinY, outputMaxY](float input) {
        return jmap(input, -24.f, 24.f, outputMinY, outputMaxY);
    };

    m_responseCurve.clear();
    m_responseCurve.startNewSubPath(responseArea.getX(), map(mags[0]));

    for (int i = 1; i < w; ++i) {
        m_responseCurve.lineTo(responseArea.getX() + i, map(mags[i]));
    }
}

// Settings positions into AudioProcessor filter
void EQGraphicComponent::mouseDrag(const juce::MouseEvent& event)
{
    int filterIndex = m_filterInterface.inBoundsOfCircle(event.x, event.y);
    if (filterIndex != -1) {
        auto& comp = m_filterInterface.fcomps[filterIndex];
        float w = comp.getWidth();
        float h = comp.getHeight();
        float newX = event.x - (w / 2);
        float newY = event.y - (h / 2);
        comp.setBounds(newX, newY, w, h);

        updateFilterParamsFromCoords(filterIndex, newX, newY);

        repaint();
    }
}

std::vector<float> EQGraphicComponent::getXs(const std::vector<float>& freqs, float left, float width)
{
    std::vector<float> xs;
    for (auto f : freqs)
    {
        auto normX = juce::mapFromLog10(f, 20.f, 20000.f);
        xs.push_back(left + width * normX);
    }

    return xs;
}

juce::Rectangle<int> EQGraphicComponent::getRenderArea()
{
    auto bounds = getLocalBounds();

    bounds.removeFromTop(12);
    bounds.removeFromBottom(2);
    bounds.removeFromLeft(20);
    bounds.removeFromRight(20);

    return bounds;
}

juce::Rectangle<int> EQGraphicComponent::getAnalysisArea()
{
    auto bounds = getRenderArea();
    bounds.removeFromTop(4);
    bounds.removeFromBottom(4);
    return bounds;
}

std::vector<float> EQGraphicComponent::getFrequencies()
{
    return std::vector<float>
    {
        20, /*30, 40,*/ 50, 100,
            200, /*300, 400,*/ 500, 1000,
            2000, /*3000, 4000,*/ 5000, 10000,
            20000
    };
}

std::vector<float> EQGraphicComponent::getGains()
{
    return std::vector<float>
    {
        -24, -12, 0, 12, 24
    };
}

void EQGraphicComponent::drawBackgroundGrid(juce::Graphics& g)
{
    using namespace juce;
    auto freqs = getFrequencies();

    auto renderArea = getAnalysisArea();

    auto left = renderArea.getX();
    auto right = renderArea.getRight();
    auto top = renderArea.getY();
    auto bottom = renderArea.getBottom();
    auto width = renderArea.getWidth();

    auto xs = getXs(freqs, left, width);

    g.setColour(Colours::dimgrey);
    for (auto x : xs)
    {
        g.drawVerticalLine(x, top, bottom);
    }

    auto gain = getGains();

    for (auto gDb : gain)
    {
        auto y = jmap(gDb, -24.f, 24.f, float(bottom), float(top));

        g.setColour(gDb == 0.f ? Colour(0u, 172u, 1u) : Colours::darkgrey);
        g.drawHorizontalLine(y, left, right);
    }
}

FilterCircle::FilterCircle()
{
    setInterceptsMouseClicks(false, false);
}

FilterCircle::~FilterCircle()
{
}

void FilterCircle::paint(juce::Graphics& g)
{
    using namespace juce;
    g.setColour(Colours::orange);

    auto bounds = getLocalBounds();
    Rectangle<float> floatBounds(bounds.getX(), bounds.getY(), bounds.getWidth()-1, bounds.getHeight()-1);
    g.drawEllipse(floatBounds, 1.0);
}
