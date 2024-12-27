#pragma once

#include <JuceHeader.h>

#include <utility>
#include <array>
#include <limits>
#include <complex>
#include <unordered_map>

/*
	Responible for:
	1) processingSample (send to biquad to process)
	2) configuring biquad with coefficients
		2b) choose coefficients based on algorithm
	3) hold state of AudioFilterParameters (cutoff, q, boost/cut)
*/

namespace dsp {

template <typename Type>
using Complex = std::complex<Type>;

inline bool checkFloatUnderflow(float& value)
{
	bool retValue = false;
	if (value > 0.0 && value < std::numeric_limits<float>::min())
	{
		value = 0;
		retValue = true;
	}
	else if (value < 0.0 && value > -std::numeric_limits<float>::min())
	{
		value = 0;
		retValue = true;
	}
	return retValue;
}

enum FilterAlgorithm {
	kLPF1P, kLPF1, kLPF2, kHPF1, kHPF2, kBPF2, kBSF2, kButterLPF2, kButterHPF2, kButterBPF2,
	kButterBSF2, kMMALPF2, kMMALPF2B, kLowShelf, kHiShelf, kNCQParaEQ, kCQParaEQ,
	kLWRLPF2, kLWRHPF2, kAPF1, kAPF2, kMatchLP2A, kMatchLP2B, kMatchBP2A,
	kMatchBP2B
};

static const std::unordered_map<FilterAlgorithm, std::string> stringToFilterAlgorithm = {
	{FilterAlgorithm::kLPF1P, "LPF-1 Only Poles"}, {FilterAlgorithm::kLPF1, "LPF-1"}, {FilterAlgorithm::kLPF2, "LPF-2"},
	{FilterAlgorithm::kHPF1, "HPF-1"}, {FilterAlgorithm::kHPF2, "HPF-2"}, {FilterAlgorithm::kBPF2, "BPF-2"},
	{FilterAlgorithm::kBSF2, "BSF-2"}, {FilterAlgorithm::kButterLPF2, "Butter LPF-2"}, {FilterAlgorithm::kButterHPF2, "Butter HPF-2"},
	{FilterAlgorithm::kButterBPF2, "Butter BPF-2"}, {FilterAlgorithm::kButterBSF2, "Butter BSF-2"}, {FilterAlgorithm::kMMALPF2, "MMA LPF-2"},
	{FilterAlgorithm::kMMALPF2B, "MMA LPF-2B"}, {FilterAlgorithm::kLowShelf, "Low Shelf"}, {FilterAlgorithm::kHiShelf, "High Shelf"},
	{FilterAlgorithm::kNCQParaEQ, "Non Const Parametric EQ"}, {FilterAlgorithm::kCQParaEQ, "Const Parametric EQ"},
	{FilterAlgorithm::kLWRLPF2, "Linkwitz Riley LPF"}, {FilterAlgorithm::kLWRHPF2, "Linkwitz Riley HPF"}, {FilterAlgorithm::kAPF1, "APF-1"},
	{FilterAlgorithm::kAPF2, "APF-2"}, {FilterAlgorithm::kMatchLP2A, "Match LPF Tight"}, {FilterAlgorithm::kMatchLP2B, "Match LPF Loose"},
	{FilterAlgorithm::kMatchBP2A, "Match BPF Tight"}, {FilterAlgorithm::kMatchBP2B, "Match BPF Loose"}
};

static juce::StringArray getFilterAlgoStrs() {
	juce::StringArray out;
	for (auto& [fa, fa_str] : stringToFilterAlgorithm) {
		out.add(fa_str);
	}
	return out;
};

class BiQuad {
public:
	BiQuad();
	~BiQuad();

	void reset();

	std::pair<float, float> getDryWet() {
		return { m_coeffs[5], m_coeffs[6] };
	}
	void setCoeffs(std::array<float, 7>& coeffs) {
		m_coeffs = std::move(coeffs);
	}
	float processSample(float xn);

	friend class IRRFilter;

private:
	std::array<float, 7> m_coeffs;
	std::array<float, 4> m_state;
};

struct FilterParams {
	float cutoffFreq;
	float q;
	float boostCutDB;
	FilterAlgorithm fa;

	bool operator==(const FilterParams& r) const {
		return std::tie(cutoffFreq, q, boostCutDB, fa) == std::tie(r.cutoffFreq, r.q, r.boostCutDB, r.fa);
	}
	bool operator!=(const FilterParams& r) const {
		return !(*this == r);
	}
};

class IRRFilter
{
public:
	IRRFilter();
	~IRRFilter();
	float processSample(float xn);
	float getMagnitudeForFrequency(float freq, float sampleRate);
	bool setCoeffs(const FilterParams& fp, float sampleRate);

private:	
	size_t getFilterOrder(float sampleRate);

	BiQuad m_biquad;
};

}