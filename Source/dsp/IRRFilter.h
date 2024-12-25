#pragma once

#include <utility>
#include <array>
#include <limits>

/*
	Responible for:
	1) processingSample (send to biquad to process)
	2) configuring biquad with coefficients
		2b) choose coefficients based on algorithm
	3) hold state of AudioFilterParameters (cutoff, q, boost/cut)
*/

namespace dsp {

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
	kLWRLPF2, kLWRHPF2, kAPF1, kAPF2, kResonA, kResonB, kMatchLP2A, kMatchLP2B, kMatchBP2A,
	kMatchBP2B, kImpInvLP1, kImpInvLP2
};

class BiQuad {
public:
	BiQuad();
	~BiQuad();

	std::pair<float, float> getDryWet() {
		return { m_coeffs[5], m_coeffs[6] };
	}
	void setCoeffs(std::array<float, 7>& coeffs) {
		m_coeffs = std::move(coeffs);
	}
	float processSample(float xn);

private:
	std::array<float, 7> m_coeffs;
	std::array<float, 4> m_state;
};

struct FilterParams {
	float cutoffFreq;
	float q;
	float boostCutDB;
};

class IRRFilter
{
public:
	IRRFilter();
	~IRRFilter();
	float processSample(float xn);

	void setFilterAlgorithm(FilterAlgorithm fa) {
		m_fa = fa;
	}
	void setSampleRate(float sampleRate) {
		m_sampleRate = sampleRate;
	}
	void setFilterParams(const FilterParams& fp) {
		m_fp = fp;
	}
	FilterAlgorithm getFilterAlgorithm() {
		return m_fa;
	}
	float getSampleRate() {
		return m_sampleRate;
	}
	FilterParams getFilterParams() {
		return m_fp;
	}

private:
	bool setCoeffs(); // Using filter params

	BiQuad m_biquad;
	FilterParams m_fp;
	float m_sampleRate;
	FilterAlgorithm m_fa;
};

}