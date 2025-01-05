/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"
#include "FilterSolver.h"

#include <algorithm>

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
    g.fillAll(Colours::grey);
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
    m_audioProcessor(ap),
    m_drawing(false)
{
    for (auto& comp : m_filterInterface.fcomps) {
        addAndMakeVisible(comp);
    }
    updateChain();
    repaint();
}

EQGraphicComponent::~EQGraphicComponent()
{
}

void EQGraphicComponent::drawTextLabels(juce::Graphics& g)
{
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
        juce::String str;
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

        juce::Rectangle<int> r;

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

        juce::Rectangle<int> r;
        r.setSize(textWidth, fontHeight);
        r.setX(1);
        r.setCentre(r.getCentreX(), y);

        g.setColour(gDb == 0.f ? juce::Colour(0u, 172u, 1u) : Colours::lightgrey);

        g.drawFittedText(str, r, juce::Justification::centredLeft, 1);
    }
}

void EQGraphicComponent::paint(juce::Graphics& g)
{
    using namespace juce;

    auto bounds = getLocalBounds();

    if (!m_drawing) {
        g.fillAll(Colours::black);
    }
    else {
        g.fillAll(Colours::darkgrey);
    }
    

    drawBackgroundGrid(g);

    drawTextLabels(g);

    g.setColour(Colours::white);
    g.strokePath(m_responseCurve, PathStrokeType(2.f));

    g.setColour(Colours::red);
    g.strokePath(m_drawCurve, PathStrokeType(2.f));

    if (!m_drawing) {
        updateResponseCurve();
    }
    
}

void EQGraphicComponent::resized() {
    auto bounds = getAnalysisArea();

    auto mapFPsToPositions = [&bounds](const auto& fp, int i){

        //DBG(width << " " << top << " " << bottom << " " << centerY);

        float midPointY = bounds.getCentreY();
        float bottom = bounds.getBottom();
        float top = bounds.getY();
        
        float posX = bounds.getX() + bounds.getWidth() * juce::mapFromLog10(fp.cutoffFreq, 20.f, 20000.f);
        float posY = 0.0;
        if (i == 1 || i == 2) { // PEAK FILTERS, TODO: Make more generic
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

    

    m_drawnPoints.resize(bounds.getWidth());
    std::fill(m_drawnPoints.begin(), m_drawnPoints.end(), -1);
    
    // TODO: Put in its own function
    for (int i = 0; i < m_filterInterface.size; ++i) {
        auto& fp = m_audioProcessor.getUserFilterParams(i);
        //DBG("index: " << i << " cutoff: " << fp.cutoffFreq << " q: " << fp.q << " fa: " << fp.fa << " boost/cut: " << fp.boostCutDB);
        juce::Rectangle<int> filterBound = mapFPsToPositions(fp, i);
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

    if (filterIndex == 1 || filterIndex == 2) {
        q = 1.f / 18.f;
    }
    m_audioProcessor.apvts.getParameter(filterID + "Q")->setValueNotifyingHost(q);
    m_audioProcessor.apvts.getParameter(filterID + "BoostCutDB")->setValueNotifyingHost(boostCutDB);

    //DBG("Cutoff: " << cutoffFreq << ", Q : " << q << ", Boost/Cut : " << boostCutDB);
}

void EQGraphicComponent::updateChain()
{
    auto filterParamsLowCut = m_audioProcessor.getUserFilterParams(0);
    auto filterParamsPeak1 = m_audioProcessor.getUserFilterParams(1);
    auto filterParamsPeak2 = m_audioProcessor.getUserFilterParams(2);
    auto filterParamsHiCut = m_audioProcessor.getUserFilterParams(3);

    //monoChain.setBypassed<ChainPositions::LowCut>(chainSettings.lowCutBypassed);
    //monoChain.setBypassed<ChainPositions::Peak>(chainSettings.peakBypassed);
    //monoChain.setBypassed<ChainPositions::HighCut>(chainSettings.highCutBypassed);

    auto lowCutCoefficients = makeLowCutFilter(filterParamsLowCut, m_audioProcessor.getSampleRate());

    auto peak1Coefficients = makePeakFilter(filterParamsPeak1, m_audioProcessor.getSampleRate());
    updateCoefficients(m_filterInterface.fchain.get<ChainPositions::Peak1>().coefficients, peak1Coefficients);

    auto peak2Coefficients = makePeakFilter(filterParamsPeak2, m_audioProcessor.getSampleRate());
    updateCoefficients(m_filterInterface.fchain.get<ChainPositions::Peak2>().coefficients, peak2Coefficients);

    auto highCutCoefficients = makeHighCutFilter(filterParamsHiCut, m_audioProcessor.getSampleRate());

    updateCutFilter(m_filterInterface.fchain.get<ChainPositions::LowCut>(),
        lowCutCoefficients,
        SLOPE);

    updateCutFilter(m_filterInterface.fchain.get<ChainPositions::HighCut>(),
        highCutCoefficients,
        SLOPE);
}

void EQGraphicComponent::updateResponseCurve()
{
    using namespace juce;
    auto responseArea = getAnalysisArea();

    auto w = responseArea.getWidth();

    auto& lowcut = m_filterInterface.fchain.get<ChainPositions::LowCut>();
    auto& peak1 = m_filterInterface.fchain.get<ChainPositions::Peak1>();
    auto& peak2 = m_filterInterface.fchain.get<ChainPositions::Peak2>();
    auto& highcut = m_filterInterface.fchain.get<ChainPositions::HighCut>();

    auto sampleRate = m_audioProcessor.getSampleRate();

    std::vector<double> mags;

    mags.resize(w);
    
    for (int i = 0; i < w; ++i) {
        double mag = 1.f;
        auto freq = mapToLog10(double(i) / double(w), 20.0, 20000.0);

        //if (!monoChain.isBypassed<ChainPositions::Peak>())
            mag *= peak1.coefficients->getMagnitudeForFrequency(freq, sampleRate);
            mag *= peak2.coefficients->getMagnitudeForFrequency(freq, sampleRate);

        //if (!monoChain.isBypassed<ChainPositions::LowCut>())
        //{
            if (!lowcut.isBypassed<0>())
                mag *= lowcut.get<0>().coefficients->getMagnitudeForFrequency(freq, sampleRate);
            if (!lowcut.isBypassed<1>())
                mag *= lowcut.get<1>().coefficients->getMagnitudeForFrequency(freq, sampleRate);
            if (!lowcut.isBypassed<2>())
                mag *= lowcut.get<2>().coefficients->getMagnitudeForFrequency(freq, sampleRate);
            if (!lowcut.isBypassed<3>())
                mag *= lowcut.get<3>().coefficients->getMagnitudeForFrequency(freq, sampleRate);
        //}

        //if (!monoChain.isBypassed<ChainPositions::HighCut>())
        //{
            if (!highcut.isBypassed<0>())
                mag *= highcut.get<0>().coefficients->getMagnitudeForFrequency(freq, sampleRate);
            if (!highcut.isBypassed<1>())
                mag *= highcut.get<1>().coefficients->getMagnitudeForFrequency(freq, sampleRate);
            if (!highcut.isBypassed<2>())
                mag *= highcut.get<2>().coefficients->getMagnitudeForFrequency(freq, sampleRate);
            if (!highcut.isBypassed<3>())
                mag *= highcut.get<3>().coefficients->getMagnitudeForFrequency(freq, sampleRate);
        //}

        mags[i] = Decibels::gainToDecibels(mag);
    }

    const float outputMinY = responseArea.getBottom();
    const float outputMaxY = responseArea.getY();
    auto map = [outputMinY, outputMaxY](float input) {
        return jmap(input, -24.f, 24.f, outputMinY, outputMaxY);
    };

    m_responseCurve.clear();
    m_responseCurve.startNewSubPath(responseArea.getX(), map(mags[0]));

    for (int i = 1; i < w; ++i) {
        m_responseCurve.lineTo(responseArea.getX() + i, map(mags[i]));
    }
}

void EQGraphicComponent::adjustFiltersAtClickPoint(int x, int y)
{
    int filterIndex = m_filterInterface.inBoundsOfCircle(x, y);
    if (filterIndex != -1) {
        auto& comp = m_filterInterface.fcomps[filterIndex];
        float w = comp.getWidth();
        float h = comp.getHeight();
        float newX = x - (w / 2);
        float newY = y - (h / 2);
        comp.setBounds(newX, newY, w, h);

        updateFilterParamsFromCoords(filterIndex, newX, newY);
    }
    updateChain();
    repaint();
}

int EQGraphicComponent::findPreviousValidX(int x) {
    x--;
    for (; x >= 0 && m_drawnPoints[x] == -1; --x);
    return x;
}

std::pair<int, int> EQGraphicComponent::findNearestAxisFromLine(int ax, int ay, int bx, int by) {
    float m = (float(by) - float(ay)) / (float(bx) - float(ax));
    float b = ay - m * ax;
    if (b > 0) return std::make_pair(0, int(b));
    return std::make_pair(int(-b / m), 0);
}

void EQGraphicComponent::resetCurveDraw(int x, int y) {
    startX = x;
    prevX = -1;
    drewToAxis = false;
    std::fill(m_drawnPoints.begin(), m_drawnPoints.end(), -1);
    m_drawCurve.clear();
    m_drawnPoints[x] = y;
    m_drawCurve.startNewSubPath(x, y);
}

void EQGraphicComponent::mouseDrag(const juce::MouseEvent& event)
{   
    if (m_drawing) {
        if (prevX < event.x) { // Drawing
            int distanceFromStart = event.x - startX;
            if (!drewToAxis && distanceFromStart > 10) {
                auto axisCrossingPoint = findNearestAxisFromLine(startX, m_drawnPoints[startX], event.x, event.y);
                m_drawnPoints[axisCrossingPoint.first] = axisCrossingPoint.second;
                m_drawCurve.lineTo(axisCrossingPoint.first, axisCrossingPoint.second);
                m_drawCurve.startNewSubPath(event.x, event.y);
                drewToAxis = true;
                // TODO: Fill in points here from currentPos to axis
            }
            m_drawnPoints[event.x] = event.y;
            m_drawCurve.lineTo(event.x, event.y);
        }
        else if (drewToAxis) { // Erasing (can only happen after we've drawn a line to the axis)
            for (int i = event.x; i < m_drawnPoints.size(); ++i) m_drawnPoints[i] = -1;
            int backAnX = findPreviousValidX(event.x);
            if (backAnX < 0) return;
            m_drawCurve.clear();
            m_drawCurve.startNewSubPath(backAnX, m_drawnPoints[backAnX]);
            for (int x = backAnX; x >= 0; --x) {
                if (m_drawnPoints[x] == -1) continue;
                m_drawCurve.lineTo(x, m_drawnPoints[x]);
            }
            m_drawCurve.startNewSubPath(backAnX, m_drawnPoints[backAnX]);
        }
        prevX = event.x;
    }
    else {
        adjustFiltersAtClickPoint(event.x, event.y);
    }
    repaint();
}

// TODO
std::vector<int> EQGraphicComponent::normalizedDrawnPoints(std::vector<int>& drawnPoints) {
    std::vector<int> normalizedPoints;
    std::copy(drawnPoints.begin(), drawnPoints.end(), std::back_inserter(normalizedPoints));
    
    // 1) Flip coords
    auto bounds = getAnalysisArea();
    for (auto& point : normalizedPoints) {
        if (point != -1) {
            point = bounds.getBottom() - point;
        }
    }

    // 2) Interpolate between points
    //auto rangeToInterp = std::make_pair(-1, normalizedPoints.size()-1);
    //int currentIndex = normalizedPoints.size();
    //while (currentIndex >= 0) {
    //    if (normalizedPoints[currentIndex] == -1) {
    //        rangeToInterp.first = currentIndex;
    //    }
    //    else {
    //        if (rangeToInterp.second != -1) {
    //            interp(rangeToInterp);
    //        }
    //        rangeToInterp.second = currentIndex;
    //    }
    //}

    return normalizedPoints;
}

void EQGraphicComponent::mouseDown(const juce::MouseEvent& event)
{
    if (m_drawing) {
        resetCurveDraw(event.x, event.y);
    }
    repaint();
}

void EQGraphicComponent::mouseUp(const juce::MouseEvent& event)
{
    // TODO: Draw line from end point to axis, also fill in drawnPoints!
    if (m_drawing) {
        DBG("RUNNING SOLVER");
        fsolve::FilterSolver fs(normalizedDrawnPoints(m_drawnPoints), m_audioProcessor.getSampleRate());
        fs.runSolver();
        DBG("DONE RUNNING");
    }
}

void EQGraphicComponent::mouseDoubleClick(const juce::MouseEvent& event)
{
    m_drawing = !m_drawing;
    if (m_drawing) {
        resetCurveDraw(event.x, event.y);
    }
    repaint();
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
    g.setColour(Colours::orange);

    auto bounds = getLocalBounds();
    juce::Rectangle<float> floatBounds(bounds.getX(), bounds.getY(), bounds.getWidth()-1, bounds.getHeight()-1);
    g.drawEllipse(floatBounds, 1.0);
}
