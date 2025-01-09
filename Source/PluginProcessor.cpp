/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
PicassoEQAudioProcessor::PicassoEQAudioProcessor()
#ifndef JucePlugin_PreferredChannelConfigurations
     : AudioProcessor (BusesProperties()
                     #if ! JucePlugin_IsMidiEffect
                      #if ! JucePlugin_IsSynth
                       .withInput  ("Input",  juce::AudioChannelSet::stereo(), true)
                      #endif
                       .withOutput ("Output", juce::AudioChannelSet::stereo(), true)
                     #endif
                       )
#endif
{
}

PicassoEQAudioProcessor::~PicassoEQAudioProcessor()
{
}

//==============================================================================
const juce::String PicassoEQAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool PicassoEQAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool PicassoEQAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool PicassoEQAudioProcessor::isMidiEffect() const
{
   #if JucePlugin_IsMidiEffect
    return true;
   #else
    return false;
   #endif
}

double PicassoEQAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int PicassoEQAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int PicassoEQAudioProcessor::getCurrentProgram()
{
    return 0;
}

void PicassoEQAudioProcessor::setCurrentProgram (int index)
{
}

const juce::String PicassoEQAudioProcessor::getProgramName (int index)
{
    return {};
}

void PicassoEQAudioProcessor::changeProgramName (int index, const juce::String& newName)
{
}

//==============================================================================
void PicassoEQAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    // Use this method as the place to do any pre-playback
    // initialisation that you need..

    juce::dsp::ProcessSpec spec;

    spec.maximumBlockSize = samplesPerBlock;

    spec.numChannels = 1;

    spec.sampleRate = sampleRate;

    leftChain.prepare(spec);
    rightChain.prepare(spec);

    updateFilters();

    osc.initialise([](float x) { return std::sin(x); });

    spec.numChannels = getTotalNumOutputChannels();
    osc.prepare(spec);
    osc.setFrequency(440);
}

void PicassoEQAudioProcessor::releaseResources()
{
    // When playback stops, you can use this as an opportunity to free up any
    // spare memory, etc.
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool PicassoEQAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
  #if JucePlugin_IsMidiEffect
    juce::ignoreUnused (layouts);
    return true;
  #else
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    // Some plugin hosts, such as certain GarageBand versions, will only
    // load plugins that support stereo bus layouts.
    if (layouts.getMainOutputChannelSet() != juce::AudioChannelSet::mono()
     && layouts.getMainOutputChannelSet() != juce::AudioChannelSet::stereo())
        return false;

    // This checks if the input layout matches the output layout
   #if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
   #endif

    return true;
  #endif
}
#endif

void PicassoEQAudioProcessor::processBlock (juce::AudioBuffer<float>& buffer, juce::MidiBuffer& midiMessages)
{
    juce::ScopedNoDenormals noDenormals;
    auto totalNumInputChannels  = getTotalNumInputChannels();
    auto totalNumOutputChannels = getTotalNumOutputChannels();

    // In case we have more outputs than inputs, this code clears any output
    // channels that didn't contain input data, (because these aren't
    // guaranteed to be empty - they may contain garbage).
    // This is here to avoid people getting screaming feedback
    // when they first compile a plugin, but obviously you don't need to keep
    // this code if your algorithm always overwrites all the output channels.
    for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
        buffer.clear (i, 0, buffer.getNumSamples());

    // This is the place where you'd normally do the guts of your plugin's
    // audio processing...
    // Make sure to reset the state if your inner loop is processing
    // the samples and the outer loop is handling the channels.
    // Alternatively, you can process the samples with the channels
    // interleaved by keeping the same state.

    updateFilters();

    juce::dsp::AudioBlock<float> block(buffer);

#ifdef OSC_DEBUG
    buffer.clear();

    for( int i = 0; i < buffer.getNumSamples(); ++i )
    {
        buffer.setSample(0, i, osc.processSample(0));
    }

    juce::dsp::ProcessContextReplacing<float> stereoContext(block);
    osc.process(stereoContext);
#endif

    
    auto leftBlock = block.getSingleChannelBlock(0);
    auto rightBlock = block.getSingleChannelBlock(1);

    juce::dsp::ProcessContextReplacing<float> leftContext(leftBlock);
    juce::dsp::ProcessContextReplacing<float> rightContext(rightBlock);

    leftChain.process(leftContext);
    rightChain.process(rightContext);

}

//==============================================================================
bool PicassoEQAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

juce::AudioProcessorEditor* PicassoEQAudioProcessor::createEditor()
{
    return new PicassoEQAudioProcessorEditor (*this);
    //return new juce::GenericAudioProcessorEditor(*this);
}

//==============================================================================
void PicassoEQAudioProcessor::getStateInformation (juce::MemoryBlock& destData)
{
    // You should use this method to store your parameters in the memory block.
    // You could do that either as raw data, or use the XML or ValueTree classes
    // as intermediaries to make it easy to save and load complex data.

    juce::MemoryOutputStream mos(destData, true);
    apvts.state.writeToStream(mos);
}

void PicassoEQAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    // You should use this method to restore your parameters from this memory block,
    // whose contents will have been created by the getStateInformation() call.

    auto tree = juce::ValueTree::readFromData(data, sizeInBytes);
    if (tree.isValid())
    {
        apvts.replaceState(tree);
        updateFilters();
    }
}

juce::AudioProcessorValueTreeState::ParameterLayout PicassoEQAudioProcessor::createParameterLayout()
{
    juce::AudioProcessorValueTreeState::ParameterLayout layout;

    for (int i = 0; i < NUM_FILTERS; ++i) {
        std::string filterName{ "Filter" + std::to_string(i) };
        float startingFreq = juce::mapToLog10(float(i + 1) / (NUM_FILTERS + 1), 20.f, 20000.f);
        layout.add(std::make_unique<juce::AudioParameterFloat>(filterName+"LowCut Freq",
            "LowCut Freq",
            juce::NormalisableRange<float>(20.f, 20000.f, 1.f, 1.f),
            startingFreq));
        layout.add(std::make_unique<juce::AudioParameterFloat>(filterName+"Q",
            "Q",
            juce::NormalisableRange<float>(0.1f, 18.f, 0.1f, 1.f),
            1.0f));
        layout.add(std::make_unique<juce::AudioParameterFloat>(filterName + "BoostCutDB",
            "BoostCutDB",
            juce::NormalisableRange<float>(-24.f, 24.f, 0.5f, 1.f),
            0.0f));
        //juce::StringArray filterAlgoStrs = dsp::getFilterAlgoStrs();
        //layout.add(std::make_unique<juce::AudioParameterChoice>(filterName+"Filter Algorithm", "Filter Algorithn", filterAlgoStrs, 0));
    }

    return layout;
}

FilterParams PicassoEQAudioProcessor::getUserFilterParams(int filterIndex) const
{
    // TODO: Cleaner way to store these names
    std::string filterName{ "Filter" + std::to_string(filterIndex) };
    FilterParams fp{};
    fp.cutoffFreq = apvts.getRawParameterValue(filterName+"LowCut Freq")->load();
    fp.q = apvts.getRawParameterValue(filterName + "Q")->load();
    fp.boostCutDB = apvts.getRawParameterValue(filterName + "BoostCutDB")->load();
    //fp.fa = static_cast<dsp::FilterAlgorithm>(apvts.getRawParameterValue(filterName + "Filter Algorithm")->load());

    return fp;
}

void updateCoefficients(Coefficients& old, const Coefficients& replacements)
{
    *old = *replacements;
}

Coefficients makePeakFilter(const FilterParams& filterParams, double sampleRate)
{
    return juce::dsp::IIR::Coefficients<float>::makePeakFilter(sampleRate,
        filterParams.cutoffFreq,
        filterParams.q,
        juce::Decibels::decibelsToGain(filterParams.boostCutDB));
}

void PicassoEQAudioProcessor::updateLowCutFilters(const FilterParams& filterParams)
{
    //auto cutCoefficients = makeLowCutFilter(filterParams, getSampleRate());
    //auto& leftLowCut = leftChain.get<ChainPositions::LowCut>();
    //auto& rightLowCut = rightChain.get<ChainPositions::LowCut>();

    //leftChain.setBypassed<ChainPositions::LowCut>(filterParams.lowCutBypassed);
    //rightChain.setBypassed<ChainPositions::LowCut>(filterParams.lowCutBypassed);

    //updateCutFilter(rightLowCut, cutCoefficients, SLOPE);
    //updateCutFilter(leftLowCut, cutCoefficients, SLOPE);
}

void PicassoEQAudioProcessor::updateHighCutFilters(const FilterParams& filterParams)
{
    //auto highCutCoefficients = makeHighCutFilter(filterParams, getSampleRate());

    //auto& leftHighCut = leftChain.get<ChainPositions::HighCut>();
    //auto& rightHighCut = rightChain.get<ChainPositions::HighCut>();

    //leftChain.setBypassed<ChainPositions::HighCut>(filterParams.highCutBypassed);
    //rightChain.setBypassed<ChainPositions::HighCut>(filterParams.highCutBypassed);

    //updateCutFilter(leftHighCut, highCutCoefficients, SLOPE);
    //updateCutFilter(rightHighCut, highCutCoefficients, SLOPE);
}

void PicassoEQAudioProcessor::updateFilters()
{
    auto fp0 = getUserFilterParams(0);
    updatePeakFilter<ChainPositions::LowCut>(fp0);

    auto fp1 = getUserFilterParams(1);
    updatePeakFilter<ChainPositions::Peak1>(fp1);

    auto fp2 = getUserFilterParams(2);
    updatePeakFilter<ChainPositions::Peak2>(fp2);

    auto fp3 = getUserFilterParams(3);
    updatePeakFilter<ChainPositions::HighCut>(fp3);
    //for (int i = 0; i < NUM_FILTERS; ++i) { // TODO: More easily scalable design (NUM_FILTERS???)
    //    auto filterParams = getUserFilterParams(i);

    //    if (i == 0) {
    //        updateLowCutFilters(filterParams);
    //    }
    //    else if (i == NUM_FILTERS - 1) {
    //        updateHighCutFilters(filterParams);
    //    }
    //    else { // TODO: Fix... so bad
    //        if (i == 1) {
    //            updatePeakFilter<ChainPositions::Peak1>(filterParams);
    //        }
    //        else if (i == 2) {
    //            updatePeakFilter<ChainPositions::Peak2>(filterParams);
    //        }
    //    } 
    //}
}

//==============================================================================
// This creates new instances of the plugin..
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new PicassoEQAudioProcessor();
}
