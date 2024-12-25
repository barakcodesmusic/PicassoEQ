
#include "IRRFilter.h"
#include <cmath>
#include <iostream>
#include <juce_core/juce_core.h>

namespace dsp {

BiQuad::BiQuad() :
	m_coeffs{0,0,0,0,0,0,0},
	m_state{0,0,0,0}
{
}

BiQuad::~BiQuad()
{
}

IRRFilter::IRRFilter()
{
}

IRRFilter::~IRRFilter()
{
}

float IRRFilter::processSample(float xn)
{
	setCoeffs();
	auto [wet, dry] = m_biquad.getDryWet();
	return dry * xn + wet * m_biquad.processSample(xn);
}

bool IRRFilter::setCoeffs()
{
	if (m_fa == FilterAlgorithm::kLPF1) {

		float theta = 2.0 * juce::MathConstants<float>::pi * m_fp.cutoffFreq / m_sampleRate;
		float gamma = std::cos(theta) / (1 + std::sin(theta));

		float a0 = (1 - gamma) / 2.0;
		float a1 = (1 - gamma) / 2.0;
		float a2 = 0.0;
		float b1 = -gamma;
		float b2 = 0.0;
		float c0 = 1.0;
		float d0 = 0.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		return true;
	} else if (m_fa == FilterAlgorithm::kButterLPF2)
	{
		// --- see book for formulae
		double theta_c = juce::MathConstants<float>::pi * m_fp.cutoffFreq / m_sampleRate;
		double C = 1.0 / std::tan(theta_c);

		// --- update coeffs
		float a0 = 1.0 / (1.0 + std::sqrt(2) * C + C * C);
		float a1 = 2.0 * a0;
		float a2 = a0;
		float b1 = 2.0 * a0 * (1.0 - C * C);
		float b2 = a0 * (1.0 - std::sqrt(2) * C + C * C);
		float c0 = 1.0;
		float d0 = 0.0;

		// --- update on calculator
		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		// --- we updated
		return true;
	}
	else {
		std::cout << "Unsupported filter algorithm: " << m_fa << std::endl;
	}
	return false;
}

float BiQuad::processSample(float xn)
{
	jassert(m_coeffs.size() == 7);
	const auto [a0, a1, a2, b1, b2, _w, _d] = m_coeffs;
	auto& [ad1, ad2, bd1, bd2] = m_state;

	float yn = a0 * xn + ad1 * a1 + ad2 * a2 - bd1 * b1 - bd2 * b2;

	checkFloatUnderflow(yn);

	ad2 = ad1;
	ad1 = xn;
	bd2 = bd1;
	bd1 = yn;

	return yn;
}

}