/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginEditor.h"
#include "PluginAgent.h"

#include <algorithm>

namespace utils {

EQPoint findNearestAxisPointRightFromLine(const juce::Rectangle<int>& bounds, const EQPoint& lineStart, const EQPoint& lineEnd) {
    float m = float(lineEnd.getY() - lineStart.getY()) / float(lineEnd.getX() - lineStart.getX());
    float b = lineStart.getY() - m * lineStart.getX();

    int axisX = bounds.getX() + bounds.getWidth();
    int axisY = int(m * axisX + b);
    if (axisY < 0) { // Cross y-axis above screen
        axisY = bounds.getCentreY(); // Default to 0db
    }
    else if (axisY >= bounds.getBottom()) { // Cross y-axis below screen
        axisX = int((bounds.getBottom() - b) / m);
        axisY = bounds.getBottom();
    }

    return { axisX, axisY };
}

EQPoint findNearestAxisPointLeftFromLine(const juce::Rectangle<int>& bounds, const EQPoint& lineStart, const EQPoint& lineEnd) {
    float m = float(lineEnd.getY() - lineStart.getY()) / float(lineEnd.getX() - lineStart.getX());
    float b = lineStart.getY() - m * lineStart.getX();

    EQPoint axis = {
        0,
        int(b)
    };
    if (b < 0) { // Cross y-axis above screen
        axis.setY(bounds.getCentreY()); // Default to 0db
    }
    else if (b >= bounds.getBottom()) { // Cross y-axis below screen
        axis = { int((bounds.getBottom() - b) / m), bounds.getBottom() };
    }
    return axis;
}

EQPoint findPreviousValidPoint(const EQPoint& startPos, const std::vector<int>& vals) {
    int x = startPos.getX() - 20; // Get a good distance
    for (; x >= 0 && (vals[x] == -1); --x);
    if (x == -1) return { -1, -1 };
    return { x, vals[x] };
}

template <typename Number>
void linearInterpolatePointToPoint(std::vector<Number>& vals, const EQPoint& start, const EQPoint& end) {

	if (start.getX() < 0 || end.getX() > vals.size()) {
		return;
	}

    jassert(end.getX() < vals.size());
    jassert(end.getX() >= start.getX());

    float xRange = end.getX() - start.x;
    float yRange = end.getY() - start.y;

    for (int i = 0; i <= xRange; ++i) {
        // Asserts protect this array access
        vals[start.x + i] = start.y + (i / xRange * yRange);
    }
}

bool validPosition(const juce::Rectangle<int>& bounds, const EQPoint& pos) {
    int l = bounds.getX();
    int r = bounds.getX() + bounds.getWidth();
    int u = bounds.getY();
    int b = bounds.getBottom();

    return pos.x >= l && pos.x <= r && pos.y >= u && pos.y <= b;
}

std::vector<float> normalizedDrawnPoints(const juce::Rectangle<int>& bounds, std::vector<int>& drawnPoints, FilterRange normalizeRange) {
    std::vector<float> normalizedPoints;
    std::copy(drawnPoints.begin(), drawnPoints.end(), std::back_inserter(normalizedPoints));

    // 1) Flip coords
    float prevNonNeg = -1;
    // TODO: Linearly interpolate between points?
    for (auto& point : normalizedPoints) {
        if (point != -1) {
            prevNonNeg = point;
        }
        else {
            point = prevNonNeg;
        }
    }
    auto it = std::find_if(normalizedPoints.begin(), normalizedPoints.end(), [](float num) {return num != -1;});
    auto firstNonNeg = *it;
    std::fill(normalizedPoints.begin(), it, firstNonNeg);

    // 2) Map to DB
    for (auto& p : normalizedPoints) {
        p = juce::jmap(p, static_cast<float>(bounds.getBottom()), static_cast<float>(bounds.getY()), GAIN_RANGE.first, GAIN_RANGE.second);
    }

    // TODO: Clean up eventually
    for (int i = 0; i < normalizeRange.startPixel; ++i) {
        normalizedPoints[i] = 0;
    }
    for (int i = normalizeRange.endPixel+1; i < normalizedPoints.size(); ++i) {
        normalizedPoints[i] = 0;
    }

    // Smooth out on left and right sides of filter range
	const int interpolateRange = 30;
    const EQPoint rightStart{ normalizeRange.endPixel, int(normalizedPoints[normalizeRange.endPixel]) };
    const EQPoint rightEnd{ normalizeRange.endPixel + interpolateRange, 0 };
    const EQPoint leftStart{ normalizeRange.startPixel - interpolateRange, 0 };
    const EQPoint leftEnd{ normalizeRange.startPixel, int(normalizedPoints[normalizeRange.startPixel]) };
    
    linearInterpolatePointToPoint(
        normalizedPoints, rightStart, rightEnd);
    linearInterpolatePointToPoint(
        normalizedPoints, leftStart, leftEnd);

    return normalizedPoints;
}

std::vector<FilterRange> calculateFilterRanges(const juce::Rectangle<int>& bounds, const std::vector<int>& drawnPoints) {

    using namespace juce;

	std::vector<FilterRange> filterRanges;

    std::vector<float> nps;
    std::copy(drawnPoints.begin(), drawnPoints.end(), std::back_inserter(nps));

    // 1) Normalize vector
    float prevNonNeg = -1;
    for (auto& point : nps) {
        if (point != -1) {
            prevNonNeg = point;
        }
        else {
            point = prevNonNeg;
        }
    }
    auto it = std::find_if(nps.begin(), nps.end(), [](float num) {return num != -1;});
    const auto firstNonNeg = *it;
    std::fill(nps.begin(), it, firstNonNeg);

    // 2) Map to DB
    for (auto p : nps) {
        p = juce::jmap(p, static_cast<float>(bounds.getBottom()), static_cast<float>(bounds.getY()), GAIN_RANGE.first, GAIN_RANGE.second);
    }

    // 3) Get direction changes for filter creation
    std::vector<int> peaksAndValleys{ 0 };
	const int distance = 20; // sensitivity to direction change
    bool directionUp = nps[distance] > nps[0];
    for (int i = distance; i < nps.size(); ++i) {
        const bool changedDirection =
            (directionUp && nps[i] < nps[i - distance]) ||
            (!directionUp && nps[i] > nps[i - distance]);
  	    if (changedDirection) {
  		    peaksAndValleys.push_back(i);
  		    directionUp = !directionUp;
  	    }
    }
    peaksAndValleys.push_back(bounds.getWidth());

    for (int i = 0; i < peaksAndValleys.size() - 1; ++i) {
		filterRanges.push_back({ peaksAndValleys[i], peaksAndValleys[i + 1] });
    }
    
    return filterRanges;
}

}

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

void EQGraphicComponent::updateDrawnCurve() {
    m_drawCurve.clear();
    int x = 0;
    for (; x < m_drawnPoints.size() && m_drawnPoints[x] == -1; ++x);
    m_drawCurve.startNewSubPath(x, m_drawnPoints[x]);
    for (; x < m_drawnPoints.size(); ++x) {
        if (m_drawnPoints[x] == -1) continue;
        m_drawCurve.lineTo(x, m_drawnPoints[x]);
    }
}

void EQGraphicComponent::paint(juce::Graphics& g)
{
    using namespace juce;

    if (!m_drawing) {
        g.fillAll(Colours::black);
        updateResponseCurve();
    }
    else {
        updateDrawnCurve();
        g.fillAll(Colours::darkgrey);
    }
    

    drawBackgroundGrid(g);

    drawTextLabels(g);

    g.setColour(Colours::white);
    g.strokePath(m_responseCurve, PathStrokeType(2.f));

    g.setColour(Colours::red);
    g.strokePath(m_drawCurve, PathStrokeType(2.f));
}

void EQGraphicComponent::setFilterAnchorPositions() {
    auto bounds = getAnalysisArea();

    auto mapFPsToPositions = [&bounds](const auto& fp, int i) {

        //DBG(width << " " << top << " " << bottom << " " << centerY);

        float midPointY = bounds.getCentreY();
        float bottom = bounds.getBottom();
        float top = bounds.getY();

        float posX = bounds.getX() + bounds.getWidth() * juce::mapFromLog10(fp.cutoffFreq, 20.f, 20000.f);
        float posY = 0.0;
        //if (i == 1 || i == 2) { // PEAK FILTERS, TODO: Make more generic
            // TODO: Figure out how to calculate q?
        posY = juce::jmap(fp.boostCutDB, GAIN_RANGE.first, GAIN_RANGE.second, bottom, top);
        //}
        //else { // FOR CUTOFF filters
        //    if (fp.q > 1) { // TODO: Not accurate midpoint?
        //        posY = juce::jmap(fp.q, 1.f, 18.f, midPointY, top); // 1 -> 18
        //    }
        //    else {
        //        posY = juce::jmap(fp.q, 0.1f, 1.f, bottom, midPointY); // .1 -> 1
        //    }
        //}

        //DBG("Posx: " << posX << " Posy: " << posY << " width: " << bounds.getWidth() << " height: " << bounds.getHeight() << " top: " << bounds.getY() << " bottom " << bounds.getBottom());

        return juce::Rectangle<int>{static_cast<int>(posX), static_cast<int>(posY), 30, 30};
        };



    m_drawnPoints.resize(getLocalBounds().getWidth());
    std::fill(m_drawnPoints.begin(), m_drawnPoints.end(), -1);

    // TODO: Put in its own function
    for (int i = 0; i < m_filterInterface.size; ++i) {
        auto& fp = m_audioProcessor.getUserFilterParams(i);
        //DBG("index: " << i << " cutoff: " << fp.cutoffFreq << " q: " << fp.q << " fa: " << fp.fa << " boost/cut: " << fp.boostCutDB);
        juce::Rectangle<int> filterBound = mapFPsToPositions(fp, i);
        m_filterInterface.fcomps[i].setBounds(filterBound);
    }

    updateChain();
}

void EQGraphicComponent::resized() {
    setFilterAnchorPositions();
}

void EQGraphicComponent::updateFilterParamsFromCoords(int filterIndex, const EQPoint& filterPos) {
    auto bounds = getAnalysisArea();

    float logMapFrequency = juce::mapToLog10(
        float(filterPos.getX()) / bounds.getWidth(),
        FREQ_RANGE.first, FREQ_RANGE.second);
    float cutoffFreq = juce::jmap(logMapFrequency, FREQ_RANGE.first, FREQ_RANGE.second, 0.f, 1.f);
    float boostCutDB = juce::jmap(
        float(filterPos.getY()),
        float(bounds.getBottom()),
        float(bounds.getY()),
        0.f,
        1.f);
    float q = 1.f / 18.f;

    DBG("Cutoff Freq: " << cutoffFreq);
    
    //if (filterPos.getY() < bounds.getCentreY()) { // In upper half
    //    q = juce::jmap(filterPos.getY(), bounds.getCentreY(), bounds.getY(), 10, 180); // 1 -> 18
    //}
    //else {
    //    q = juce::jmap(filterPos.getX(), bounds.getBottom(), bounds.getCentreY(), 1, 10); // .1 -> 1
    //}

    m_audioProcessor.apvts.getParameter(cutoffParamFromIndex(filterIndex))->setValueNotifyingHost(cutoffFreq);
    m_audioProcessor.apvts.getParameter(gainDBParamFromIndex(filterIndex))->setValueNotifyingHost(boostCutDB);
    m_audioProcessor.apvts.getParameter(qParamFromIndex(filterIndex))->setValueNotifyingHost(q);
}

void EQGraphicComponent::updateChain()
{
    //auto filterParamsLowCut = m_audioProcessor.getUserFilterParams(0);
    auto filterParamsPeak0 = m_audioProcessor.getUserFilterParams(0);
    auto filterParamsPeak1 = m_audioProcessor.getUserFilterParams(1);
    auto filterParamsPeak2 = m_audioProcessor.getUserFilterParams(2);
    auto filterParamsPeak3 = m_audioProcessor.getUserFilterParams(3);
    //auto filterParamsHiCut = m_audioProcessor.getUserFilterParams(3);

    //monoChain.setBypassed<ChainPositions::LowCut>(chainSettings.lowCutBypassed);
    //monoChain.setBypassed<ChainPositions::Peak>(chainSettings.peakBypassed);
    //monoChain.setBypassed<ChainPositions::HighCut>(chainSettings.highCutBypassed);

    //auto lowCutCoefficients = makeLowCutFilter(filterParamsLowCut, m_audioProcessor.getSampleRate());

    auto peak0Coefficients = makePeakFilter(filterParamsPeak0, m_audioProcessor.getSampleRate());
    updateCoefficients(m_filterInterface.fchain.get<0>().coefficients, peak0Coefficients);

    auto peak1Coefficients = makePeakFilter(filterParamsPeak1, m_audioProcessor.getSampleRate());
    updateCoefficients(m_filterInterface.fchain.get<ChainPositions::Peak1>().coefficients, peak1Coefficients);

    auto peak2Coefficients = makePeakFilter(filterParamsPeak2, m_audioProcessor.getSampleRate());
    updateCoefficients(m_filterInterface.fchain.get<ChainPositions::Peak2>().coefficients, peak2Coefficients);

    auto peak3Coefficients = makePeakFilter(filterParamsPeak3, m_audioProcessor.getSampleRate());
    updateCoefficients(m_filterInterface.fchain.get<3>().coefficients, peak3Coefficients);

    //auto highCutCoefficients = makeHighCutFilter(filterParamsHiCut, m_audioProcessor.getSampleRate());

    //updateCutFilter(m_filterInterface.fchain.get<ChainPositions::LowCut>(),
    //    lowCutCoefficients,
    //    SLOPE);

    //updateCutFilter(m_filterInterface.fchain.get<ChainPositions::HighCut>(),
    //    highCutCoefficients,
    //    SLOPE);
}

void EQGraphicComponent::updateResponseCurve()
{
    using namespace juce;
    auto responseArea = getAnalysisArea();

    auto w = responseArea.getWidth();

    //auto& lowcut = m_filterInterface.fchain.get<ChainPositions::LowCut>();
    auto& peak0 = m_filterInterface.fchain.get<0>();
    auto& peak1 = m_filterInterface.fchain.get<ChainPositions::Peak1>();
    auto& peak2 = m_filterInterface.fchain.get<ChainPositions::Peak2>();
    auto& peak3 = m_filterInterface.fchain.get<3>();
    //auto& highcut = m_filterInterface.fchain.get<ChainPositions::HighCut>();

    auto sampleRate = m_audioProcessor.getSampleRate();

    std::vector<float> mags;

    mags.resize(w);
    
    for (int i = 0; i < w; ++i) {
        float mag = 1.f;
        float freq = mapToLog10(i / float(w), FREQ_RANGE.first, FREQ_RANGE.second);

        //if (!monoChain.isBypassed<ChainPositions::Peak>())
            mag *= peak0.coefficients->getMagnitudeForFrequency(freq, sampleRate);
            mag *= peak1.coefficients->getMagnitudeForFrequency(freq, sampleRate);
            mag *= peak2.coefficients->getMagnitudeForFrequency(freq, sampleRate);
            mag *= peak3.coefficients->getMagnitudeForFrequency(freq, sampleRate);

        //if (!monoChain.isBypassed<ChainPositions::LowCut>())
        //{
            //if (!lowcut.isBypassed<0>())
            //    mag *= lowcut.get<0>().coefficients->getMagnitudeForFrequency(freq, sampleRate);
            //if (!lowcut.isBypassed<1>())
            //    mag *= lowcut.get<1>().coefficients->getMagnitudeForFrequency(freq, sampleRate);
            //if (!lowcut.isBypassed<2>())
            //    mag *= lowcut.get<2>().coefficients->getMagnitudeForFrequency(freq, sampleRate);
            //if (!lowcut.isBypassed<3>())
            //    mag *= lowcut.get<3>().coefficients->getMagnitudeForFrequency(freq, sampleRate);
        //}

        //if (!monoChain.isBypassed<ChainPositions::HighCut>())
        //{
            //if (!highcut.isBypassed<0>())
            //    mag *= highcut.get<0>().coefficients->getMagnitudeForFrequency(freq, sampleRate);
            //if (!highcut.isBypassed<1>())
            //    mag *= highcut.get<1>().coefficients->getMagnitudeForFrequency(freq, sampleRate);
            //if (!highcut.isBypassed<2>())
            //    mag *= highcut.get<2>().coefficients->getMagnitudeForFrequency(freq, sampleRate);
            //if (!highcut.isBypassed<3>())
            //    mag *= highcut.get<3>().coefficients->getMagnitudeForFrequency(freq, sampleRate);
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

void EQGraphicComponent::adjustFiltersAtClickPoint(const EQPoint& filterPos)
{
    int filterIndex = m_filterInterface.inBoundsOfCircle(filterPos);
    if (filterIndex != -1) {
        auto& comp = m_filterInterface.fcomps[filterIndex];
        int newFilterX = filterPos.getX() - (comp.getWidth() / 2);
        int newFilterY = filterPos.getY() - (comp.getHeight() / 2);
        comp.setBounds(newFilterX, newFilterY, comp.getWidth(), comp.getHeight());

        updateFilterParamsFromCoords(filterIndex, {newFilterX - comp.getWidth()/3, newFilterY - comp.getHeight()});
    }
    updateChain();
    repaint();
}

void EQGraphicComponent::resetCurveDraw(const EQPoint& mouseDownPoint) {
    m_prevX = -1;
    m_drewToAxis = false;
    std::fill(m_drawnPoints.begin(), m_drawnPoints.end(), -1);
    m_drawnPoints[mouseDownPoint.getX()] = mouseDownPoint.getY();
    m_drawCurve.clear();
    m_drawCurve.startNewSubPath(mouseDownPoint.getX(), mouseDownPoint.getY());

    // Stop running threads
    threadManager.reset();
}

void EQGraphicComponent::mouseDrag(const juce::MouseEvent& event)
{   
    if (m_drawing) {
        // Draw view
        if (m_prevX < event.x) { // Moving mouse right is drawing

            EQPoint mousePos{ event.x, event.y };
            if (!utils::validPosition(getAnalysisArea(), mousePos)) return;

            EQPoint previousPoint = utils::findPreviousValidPoint(mousePos, m_drawnPoints);
            if (!m_drewToAxis && previousPoint.getX() != -1) {
                auto axis = utils::findNearestAxisPointLeftFromLine(getAnalysisArea(), previousPoint, mousePos);
                utils::linearInterpolatePointToPoint(m_drawnPoints, axis, mousePos);
                m_drewToAxis = true;
            }
            m_drawnPoints[event.x] = event.y;
        }
        else if (m_drewToAxis) { // Erasing (can only happen after we've drawn a line to the axis)
            for (int i = event.x; i < m_drawnPoints.size(); ++i) m_drawnPoints[i] = -1;
        }
        m_prevX = event.x;
    }
    else {
        // Filter view
        adjustFiltersAtClickPoint({event.x, event.y});
    }
    repaint();
}

void EQGraphicComponent::mouseDown(const juce::MouseEvent& event)
{
    EQPoint mouseDownPoint{ event.x, event.y };
    if (!utils::validPosition(getAnalysisArea(), mouseDownPoint)) return;

    if (m_drawing) {
        resetCurveDraw(mouseDownPoint);
    }
    repaint();
}

void EQGraphicComponent::solverThreadExitedCallback(const int filterIndex, const FilterParams fp) {
    // TODO: Pass down index, then derive name from there
    m_audioProcessor.apvts.getParameter(cutoffParamFromIndex(filterIndex))->setValueNotifyingHost(mapFreqToFrac(fp.cutoffFreq));
    m_audioProcessor.apvts.getParameter(qParamFromIndex(filterIndex))->setValueNotifyingHost(mapQToFrac(fp.q));
    m_audioProcessor.apvts.getParameter(gainDBParamFromIndex(filterIndex))->setValueNotifyingHost(mapDBToFrac(fp.boostCutDB));

    setFilterAnchorPositions();
    repaint();
}

void EQGraphicComponent::mouseUp(const juce::MouseEvent& event)
{
    EQPoint mouseUpPoint{ event.x, event.y };
    if (!utils::validPosition(getAnalysisArea(), mouseUpPoint)) return;

    if (m_drawing) { 
        EQPoint previousPoint = utils::findPreviousValidPoint({event.x, event.y}, m_drawnPoints);
        if (previousPoint.getX() >= 0) {
            const auto axis = utils::findNearestAxisPointRightFromLine(getAnalysisArea(), previousPoint, mouseUpPoint);

            utils::linearInterpolatePointToPoint(m_drawnPoints, mouseUpPoint, axis);
        }

        auto bounds = getAnalysisArea();
		std::vector<FilterRange> filterRanges = utils::calculateFilterRanges(bounds, m_drawnPoints);

        auto cb = std::bind(&EQGraphicComponent::solverThreadExitedCallback, this, std::placeholders::_1, std::placeholders::_2);
        
        // TODO: Need to keep track of filter positions and create filter chain based on this...
        for (int i = 0; i < filterRanges.size(); ++i) {
            const auto& filterRange = filterRanges[i];
            threadManager.addThread(
                i,
                utils::normalizedDrawnPoints(bounds, m_drawnPoints, filterRange),
                makePeakFilter, // TODO: figure out appropriate filter
                filterRange.getFilterParamsGuess(), // TODO: better starting conditions
                m_audioProcessor.getSampleRate(),
                cb
            );
        }
   
        for (const auto& thread : threadManager.getThreads()) {
	 	    thread->startThread();
        }

        m_drawing = false;
    }
    repaint();
}

void EQGraphicComponent::mouseDoubleClick(const juce::MouseEvent& event)
{
    pagent::PluginAgent agent(pagent::PluginType::Compressor);
	std::string inference = agent.inferParameters(std::string("test"));
	DBG(inference);
    //m_drawing = !m_drawing;
    //if (m_drawing) {
    //    resetCurveDraw({event.x, event.y});
    //}
    //repaint();
}

std::vector<float> EQGraphicComponent::getXs(const std::vector<float>& freqs, float left, float width)
{
    std::vector<float> xs;
    for (auto f : freqs)
    {
        auto normX = juce::mapFromLog10(f, FREQ_RANGE.first, FREQ_RANGE.second);
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

    g.setColour(Colours::dimgrey);

    auto renderArea = getAnalysisArea();
    auto xs = getXs(getFrequencies(), renderArea.getX(), renderArea.getWidth());
    for (auto x : xs)
    {
        g.drawVerticalLine(x, renderArea.getY(), renderArea.getBottom());
    }

    for (auto gDb : getGains())
    {
        auto y = jmap(gDb, -24.f, 24.f, float(renderArea.getBottom()), float(renderArea.getY()));

        g.setColour(gDb == 0.f ? Colour(0u, 172u, 1u) : Colours::darkgrey);
        g.drawHorizontalLine(y, renderArea.getX(), renderArea.getRight());
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
    juce::Rectangle<float> floatBounds(bounds.getX(), bounds.getY(), bounds.getWidth()-5, bounds.getHeight()-5);
    g.drawEllipse(floatBounds, 1.0);
}
