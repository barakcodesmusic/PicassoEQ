
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
	} else if (m_fa == FilterAlgorithm::kHPF1)
	{
		float theta = 2.0 * juce::MathConstants<float>::pi * m_fp.cutoffFreq / m_sampleRate;
		float gamma = std::cos(theta) / (1 + std::sin(theta));

		float a0 = (1 + gamma) / 2.0;
		float a1 = (1 + gamma) / 2.0;
		float a2 = 0.0;
		float b1 = -gamma;
		float b2 = 0.0;
		float c0 = 1.0;
		float d0 = 0.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		return true;
	} else if (m_fa == FilterAlgorithm::kLPF2)
	{
		float theta = 2.0 * juce::MathConstants<float>::pi * m_fp.cutoffFreq / m_sampleRate;
		float d = 1 / m_fp.q;

		float beta = 0.5 * ((1 - (d / 2) * std::sin(theta)) / (1 + (d / 2) * std::sin(theta)));
		float gamma = (0.5 + beta) * std::cos(theta);

		float a0 = (0.5 + beta - gamma) / 2.0;
		float a1 = 0.5 + beta - gamma;
		float a2 = (0.5 + beta - gamma) / 2.0;
		float b1 = -2*gamma;
		float b2 = 2*beta;
		float c0 = 1.0;
		float d0 = 0.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		return true;
	}
	else if (m_fa == FilterAlgorithm::kHPF2)
	{
		float theta = 2.0 * juce::MathConstants<float>::pi * m_fp.cutoffFreq / m_sampleRate;
		float d = 1 / m_fp.q;

		float beta = 0.5 * ((1 - (d / 2) * std::sin(theta)) / (1 + (d / 2) * std::sin(theta)));
		float gamma = (0.5 + beta) * std::cos(theta);

		float a0 = (0.5 + beta + gamma) / 2.0;
		float a1 = -(0.5 + beta + gamma);
		float a2 = (0.5 + beta + gamma) / 2.0;
		float b1 = -2 * gamma;
		float b2 = 2 * beta;
		float c0 = 1.0;
		float d0 = 0.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		return true;
	} else if (m_fa == FilterAlgorithm::kBPF2) {
		float k = std::tan(juce::MathConstants<float>::pi * m_fp.cutoffFreq / m_sampleRate);
		float delta = std::pow(k, 2) * m_fp.q + k + m_fp.q;
		float a0 = k / delta;
		float a1 = 0.0;
		float a2 = -k / delta;
		float b1 = 2 * m_fp.q * (std::pow(k, 2) - 1) / delta;
		float b2 = std::pow(k, 2) * m_fp.q - k + m_fp.q / delta;
		float c0 = 1.0;
		float d0 = 0.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		return true;
	} else if (m_fa == FilterAlgorithm::kBSF2) {
		float k = std::tan(juce::MathConstants<float>::pi * m_fp.cutoffFreq / m_sampleRate);
		float delta = std::pow(k, 2) * m_fp.q + k + m_fp.q;
		float a0 = m_fp.q * (std::pow(k,2)+1) / delta;
		float a1 = 2 * m_fp.q * (std::pow(k, 2) - 1) / delta;
		float a2 = m_fp.q * (std::pow(k,2)+1) / delta;
		float b1 = 2 * m_fp.q * (std::pow(k, 2) - 1) / delta;
		float b2 = std::pow(k, 2) * m_fp.q - k + m_fp.q / delta;
		float c0 = 1.0;
		float d0 = 0.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		return true;
	} else if (m_fa == FilterAlgorithm::kButterLPF2) {
		float theta_c = juce::MathConstants<float>::pi * m_fp.cutoffFreq / m_sampleRate;
		float C = 1.0 / std::tan(theta_c);

		float a0 = 1.0 / (1.0 + std::sqrt(2) * C + std::pow(C,2));
		float a1 = 2.0 * a0;
		float a2 = a0;
		float b1 = 2.0 * a0 * (1.0 - std::pow(C, 2));
		float b2 = a0 * (1.0 - std::sqrt(2) * C + std::pow(C, 2));
		float c0 = 1.0;
		float d0 = 0.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		return true;
	}
	else if (m_fa == FilterAlgorithm::kButterHPF2) {
		float theta_c = juce::MathConstants<float>::pi * m_fp.cutoffFreq / m_sampleRate;
		float C = std::tan(theta_c);

		float a0 = 1.0 / (1.0 + std::sqrt(2) * C + std::pow(C, 2));
		float a1 = -2.0 * a0;
		float a2 = a0;
		float b1 = 2.0 * a0 * (std::pow(C, 2) - 1.0);
		float b2 = a0 * (1.0 - std::sqrt(2) * C + std::pow(C, 2));
		float c0 = 1.0;
		float d0 = 0.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		return true;
	}
	else if (m_fa == FilterAlgorithm::kButterBPF2) {
		float BW = m_fp.cutoffFreq / m_fp.q;
		float C = 1/std::tan(juce::MathConstants<float>::pi * m_fp.cutoffFreq * BW / m_sampleRate);
		float D = 2 * std::cos(2 * juce::MathConstants<float>::pi * m_fp.cutoffFreq / m_sampleRate);

		float a0 = 1.0 / (1.0 + C);
		float a1 = 0.0;
		float a2 = -a0;
		float b1 = -a0 * C * D;
		float b2 = a0 * (C-1);
		float c0 = 1.0;
		float d0 = 0.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		return true;
	}
	else if (m_fa == FilterAlgorithm::kButterBSF2) {
		float BW = m_fp.cutoffFreq / m_fp.q;
		float C = std::tan(juce::MathConstants<float>::pi * m_fp.cutoffFreq * BW / m_sampleRate);
		float D = 2 * std::cos(2 * juce::MathConstants<float>::pi * m_fp.cutoffFreq / m_sampleRate);

		float a0 = 1.0 / (1.0 + C);
		float a1 = -a0*D;
		float a2 = a0;
		float b1 = -a0 * D;
		float b2 = a0 * (1 - C);
		float c0 = 1.0;
		float d0 = 0.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		return true;
	}
	else if (m_fa == FilterAlgorithm::kLWRLPF2) {
		float theta_c = juce::MathConstants<float>::pi * m_fp.cutoffFreq / m_sampleRate;
		float omega_c = juce::MathConstants<float>::pi * m_fp.cutoffFreq;
		float kappa = omega_c / std::tan(theta_c);
		float delta = std::pow(kappa, 2) + std::pow(omega_c, 2) + 2 * kappa * omega_c;

		float a0 = std::pow(omega_c, 2) / delta;
		float a1 = 2.0 * std::pow(omega_c, 2) / delta;
		float a2 = std::pow(omega_c, 2) / delta;
		float b1 = (- 2.0 * std::pow(kappa, 2) + 2.0 * std::pow(omega_c, 2)) / delta;
		float b2 = (-2.0 * kappa * omega_c + std::pow(kappa, 2) + std::pow(omega_c, 2)) / delta;
		float c0 = 1.0;
		float d0 = 0.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		return true;
	}
	else if (m_fa == FilterAlgorithm::kLWRHPF2) {
		float theta_c = juce::MathConstants<float>::pi * m_fp.cutoffFreq / m_sampleRate;
		float omega_c = juce::MathConstants<float>::pi * m_fp.cutoffFreq;
		float kappa = omega_c / std::tan(theta_c);
		float delta = std::pow(kappa, 2) + std::pow(omega_c, 2) + 2 * kappa * omega_c;

		float a0 = std::pow(kappa, 2) / delta;
		float a1 = -2.0 * std::pow(kappa, 2) / delta;
		float a2 = std::pow(omega_c, 2) / delta;
		float b1 = (-2.0 * std::pow(kappa, 2) + 2.0 * std::pow(omega_c, 2)) / delta;
		float b2 = (-2.0 * kappa * omega_c + std::pow(kappa, 2) + std::pow(omega_c, 2)) / delta;
		float c0 = 1.0;
		float d0 = 0.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		return true;
	}
	else if (m_fa == FilterAlgorithm::kAPF1) {
		float alpha = (std::tan(juce::MathConstants<float>::pi * m_fp.cutoffFreq / m_sampleRate) - 1) / (std::tan(juce::MathConstants<float>::pi * m_fp.cutoffFreq / m_sampleRate) + 1);
		float a0 = alpha;
		float a1 = 1.0;
		float a2 = 0.0;
		float b1 = alpha;
		float b2 = 0.0;
		float c0 = 1.0;
		float d0 = 0.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		return true;
	}
	else if (m_fa == FilterAlgorithm::kAPF2) {
		float BW = m_fp.cutoffFreq / m_fp.q;
		float alpha = (std::tan(juce::MathConstants<float>::pi * BW / m_sampleRate) - 1) / (std::tan(juce::MathConstants<float>::pi * BW / m_sampleRate) + 1);
		float beta = -std::cos(2 * juce::MathConstants<float>::pi * m_fp.cutoffFreq / m_sampleRate);
		float a0 = alpha;
		float a1 = beta * (1 - alpha);
		float a2 = 1.0;
		float b1 = beta * (1 - alpha);
		float b2 = -alpha;
		float c0 = 1.0;
		float d0 = 0.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		return true;
	}
	else if (m_fa == FilterAlgorithm::kLowShelf) {
		float theta_c = 2*juce::MathConstants<float>::pi * m_fp.cutoffFreq / m_sampleRate;
		float mu = std::pow(10, m_fp.boostCutDB / 20);
		float beta = 4 / (1 + mu);
		float delta = beta * std::tan(theta_c / 2);
		float gamma = (1 - delta) / (1 + delta);
		float a0 = (1 - gamma) / 2;
		float a1 = (1 - gamma) / 2;
		float a2 = 0.0;
		float b1 = -gamma;
		float b2 = 0.0;
		float c0 = mu - 1.0;
		float d0 = 1.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		return true;
	}
	else if (m_fa == FilterAlgorithm::kHiShelf) {
		float theta_c = 2 * juce::MathConstants<float>::pi * m_fp.cutoffFreq / m_sampleRate;
		float mu = std::pow(10, m_fp.boostCutDB / 20);
		float beta = (1 + mu) / 4;
		float delta = beta * std::tan(theta_c / 2);
		float gamma = (1 - delta) / (1 + delta);
		float a0 = (1 + gamma) / 2;
		float a1 = -(1 + gamma) / 2;
		float a2 = 0.0;
		float b1 = -gamma;
		float b2 = 0.0;
		float c0 = mu - 1.0;
		float d0 = 1.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		return true;
	}
	else if (m_fa == FilterAlgorithm::kNCQParaEQ) {
		float theta_c = 2 * juce::MathConstants<float>::pi * m_fp.cutoffFreq / m_sampleRate;
		float mu = std::pow(10, m_fp.boostCutDB / 20);
		float zeta = 4 / (1 + mu);
		float beta = 0.5 * (1-zeta*std::tan(theta_c / 2*m_fp.q))/(1 + zeta * std::tan(theta_c / 2 * m_fp.q));
		float delta = (0.5 + beta) * std::cos(theta_c);
		float a0 = 0.5 - beta;
		float a1 = 0.0;
		float a2 = -(0.5 - beta);
		float b1 = -delta;
		float b2 = 2*beta;
		float c0 = mu - 1.0;
		float d0 = 1.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		return true;
	}
	else if (m_fa == FilterAlgorithm::kCQParaEQ) {
		float K = std::tan(juce::MathConstants<float>::pi * m_fp.cutoffFreq / m_sampleRate);
		float Vo = pow(10.0, m_fp.boostCutDB / 20.0);
		bool bBoost = m_fp.boostCutDB >= 0 ? true : false;

		float d_0 = 1.0 + (1.0 / m_fp.q) * K + K * K;
		float e_0 = 1.0 + (1.0 / (Vo * m_fp.q)) * K + K * K;
		float alpha = 1.0 + (Vo / m_fp.q) * K + K * K;
		float beta = 2.0 * (K * K - 1.0);
		float gamma = 1.0 - (Vo / m_fp.q) * K + K * K;
		float delta = 1.0 - (1.0 / m_fp.q) * K + K * K;
		float eta = 1.0 - (1.0 / (Vo * m_fp.q)) * K + K * K;

		float a0 = bBoost ? alpha / d_0 : d_0 / e_0;
		float a1 = bBoost ? beta / d_0 : beta / e_0;
		float a2 = bBoost ? gamma / d_0 : delta / e_0;
		float b1 = bBoost ? beta / d_0 : beta / e_0;
		float b2 = bBoost ? delta / d_0 : eta / e_0;
		float c0 = 1.0;
		float d0 = 0.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		return true;
	}
	else if (m_fa == FilterAlgorithm::kLPF1P) {
		double theta_c = 2.0 * juce::MathConstants<float>::pi * m_fp.cutoffFreq / m_sampleRate;
		double gamma = 2.0 - cos(theta_c);

		double filter_b1 = pow((gamma * gamma - 1.0), 0.5) - gamma;
		double filter_a0 = 1.0 + filter_b1;

		float a0 = filter_a0;
		float a1 = 0.0;
		float a2 = 0.0;
		float b1 = filter_b1;
		float b2 = 0.0;
		float c0 = 1.0;
		float d0 = 0.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		return true;
	}
	else if (m_fa == FilterAlgorithm::kMMALPF2 || m_fa == FilterAlgorithm::kMMALPF2B) {
		double theta_c = 2.0 * juce::MathConstants<float>::pi * m_fp.cutoffFreq / m_sampleRate;
		double resonance_dB = 0;

		if (m_fp.q > 0.707)
		{
			double peak = m_fp.q * m_fp.q / pow(m_fp.q * m_fp.q - 0.25, 0.5);
			resonance_dB = 20.0 * log10(peak);
		}

		double resonance = (cos(theta_c) + (sin(theta_c) * sqrt(pow(10.0, (resonance_dB / 10.0)) - 1))) / ((pow(10.0, (resonance_dB / 20.0)) * sin(theta_c)) + 1);
		double g = pow(10.0, (-resonance_dB / 40.0));

		if (m_fa == FilterAlgorithm::kMMALPF2B)
			g = 1.0;

		double filter_b1 = (-2.0) * resonance * cos(theta_c);
		double filter_b2 = resonance * resonance;
		double filter_a0 = g * (1 + filter_b1 + filter_b2);

		float a0 = filter_a0;
		float a1 = 0.0;
		float a2 = 0.0;
		float b1 = filter_b1;
		float b2 = filter_b2;
		float c0 = 1.0;
		float d0 = 0.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		return true;
	}
	// --- kMatchLP2A = TIGHT fit LPF vicanek algo
	else if (m_fa == FilterAlgorithm::kMatchLP2A)
	{
		// http://vicanek.de/articles/BiquadFits.pdf
		double theta_c = 2.0 * juce::MathConstants<float>::pi * m_fp.cutoffFreq/m_sampleRate;

		double q = 1.0 / (2.0 * m_fp.q);

		// --- impulse invariant
		double b_1 = 0.0;
		double b_2 = exp(-2.0 * q * theta_c);
		if (q <= 1.0)
		{
			b_1 = -2.0 * exp(-q * theta_c) * cos(pow((1.0 - q * q), 0.5) * theta_c);
		}
		else
		{
			b_1 = -2.0 * exp(-q * theta_c) * cosh(pow((q * q - 1.0), 0.5) * theta_c);
		}

		// --- TIGHT FIT --- //
		double B0 = (1.0 + b_1 + b_2) * (1.0 + b_1 + b_2);
		double B1 = (1.0 - b_1 + b_2) * (1.0 - b_1 + b_2);
		double B2 = -4.0 * b_2;

		double phi_0 = 1.0 - sin(theta_c / 2.0) * sin(theta_c / 2.0);
		double phi_1 = sin(theta_c / 2.0) * sin(theta_c / 2.0);
		double phi_2 = 4.0 * phi_0 * phi_1;

		double R1 = (B0 * phi_0 + B1 * phi_1 + B2 * phi_2) * (m_fp.q * m_fp.q);
		double A0 = B0;
		double A1 = (R1 - A0 * phi_0) / phi_1;

		if (A0 < 0.0)
			A0 = 0.0;
		if (A1 < 0.0)
			A1 = 0.0;

		double a_0 = 0.5 * (pow(A0, 0.5) + pow(A1, 0.5));
		double a_1 = pow(A0, 0.5) - a_0;
		double a_2 = 0.0;

		float a0 = a_0;
		float a1 = a_1;
		float a2 = a_2;
		float b1 = b_1;
		float b2 = b_2;
		float c0 = 1.0;
		float d0 = 0.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		return true;
		}
		// --- kMatchLP2B = LOOSE fit LPF vicanek algo
	else if (m_fa == FilterAlgorithm::kMatchLP2B)
	{
		// http://vicanek.de/articles/BiquadFits.pdf
		double theta_c = 2.0 * juce::MathConstants<float>::pi * m_fp.cutoffFreq/m_sampleRate;
		double q = 1.0 / (2.0 * m_fp.q);

		// --- impulse invariant
		double b_1 = 0.0;
		double b_2 = exp(-2.0 * q * theta_c);
		if (q <= 1.0)
		{
			b_1 = -2.0 * exp(-q * theta_c) * cos(pow((1.0 - q * q), 0.5) * theta_c);
		}
		else
		{
			b_1 = -2.0 * exp(-q * theta_c) * cosh(pow((q * q - 1.0), 0.5) * theta_c);
		}

		// --- LOOSE FIT --- //
		double f0 = theta_c / juce::MathConstants<float>::pi; // note f0 = fraction of pi, so that f0 = 1.0 = pi = Nyquist

		double r0 = 1.0 + b_1 + b_2;
		double denom = (1.0 - f0 * f0) * (1.0 - f0 * f0) + (f0 * f0) / (m_fp.q * m_fp.q);
		denom = pow(denom, 0.5);
		double r1 = ((1.0 - b_1 + b_2) * f0 * f0) / (denom);

		double a_0 = (r0 + r1) / 2.0;
		double a_1 = r0 - a_0;
		double a_2 = 0.0;

		float a0 = a_0;
		float a1 = a_1;
		float a2 = a_2;
		float b1 = b_1;
		float b2 = b_2;
		float c0 = 1.0;
		float d0 = 0.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		return true;
		}
		// --- kMatchBP2A = TIGHT fit BPF vicanek algo
	else if (m_fa == FilterAlgorithm::kMatchBP2A)
	{
		// http://vicanek.de/articles/BiquadFits.pdf
		double theta_c = 2.0 * juce::MathConstants<float>::pi * m_fp.cutoffFreq/m_sampleRate;
		double q = 1.0 / (2.0 * m_fp.q);

		// --- impulse invariant
		double b_1 = 0.0;
		double b_2 = exp(-2.0 * q * theta_c);
		if (q <= 1.0)
		{
			b_1 = -2.0 * exp(-q * theta_c) * cos(pow((1.0 - q * q), 0.5) * theta_c);
		}
		else
		{
			b_1 = -2.0 * exp(-q * theta_c) * cosh(pow((q * q - 1.0), 0.5) * theta_c);
		}

		// --- TIGHT FIT --- //
		double B0 = (1.0 + b_1 + b_2) * (1.0 + b_1 + b_2);
		double B1 = (1.0 - b_1 + b_2) * (1.0 - b_1 + b_2);
		double B2 = -4.0 * b_2;

		double phi_0 = 1.0 - sin(theta_c / 2.0) * sin(theta_c / 2.0);
		double phi_1 = sin(theta_c / 2.0) * sin(theta_c / 2.0);
		double phi_2 = 4.0 * phi_0 * phi_1;

		double R1 = B0 * phi_0 + B1 * phi_1 + B2 * phi_2;
		double R2 = -B0 + B1 + 4.0 * (phi_0 - phi_1) * B2;

		double A2 = (R1 - R2 * phi_1) / (4.0 * phi_1 * phi_1);
		double A1 = R2 + 4.0 * (phi_1 - phi_0) * A2;

		double a_1 = -0.5 * (pow(A1, 0.5));
		double a_0 = 0.5 * (pow((A2 + (a_1 * a_1)), 0.5) - a_1);
		double a_2 = -a_0 - a_1;

		float a0 = a_0;
		float a1 = a_1;
		float a2 = a_2;
		float b1 = b_1;
		float b2 = b_2;
		float c0 = 1.0;
		float d0 = 0.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

		return true;
		}
		// --- kMatchBP2B = LOOSE fit BPF vicanek algo
	else if (m_fa == FilterAlgorithm::kMatchBP2B)
	{
		// http://vicanek.de/articles/BiquadFits.pdf
		double theta_c = 2.0 * juce::MathConstants<float>::pi * m_fp.cutoffFreq/m_sampleRate;
		double q = 1.0 / (2.0 * m_fp.q);

		// --- impulse invariant
		double b_1 = 0.0;
		double b_2 = exp(-2.0 * q * theta_c);
		if (q <= 1.0)
		{
			b_1 = -2.0 * exp(-q * theta_c) * cos(pow((1.0 - q * q), 0.5) * theta_c);
		}
		else
		{
			b_1 = -2.0 * exp(-q * theta_c) * cosh(pow((q * q - 1.0), 0.5) * theta_c);
		}

		// --- LOOSE FIT --- //
		double f0 = theta_c / juce::MathConstants<float>::pi; // note f0 = fraction of pi, so that f0 = 1.0 = pi = Nyquist

		double r0 = (1.0 + b_1 + b_2) / (juce::MathConstants<float>::pi * f0 * m_fp.q);
		double denom = (1.0 - f0 * f0) * (1.0 - f0 * f0) + (f0 * f0) / (m_fp.q * m_fp.q);
		denom = pow(denom, 0.5);

		double r1 = ((1.0 - b_1 + b_2) * (f0 / m_fp.q)) / (denom);

		double a_1 = -r1 / 2.0;
		double a_0 = (r0 - a_1) / 2.0;
		double a_2 = -a_0 - a_1;

		float a0 = a_0;
		float a1 = a_1;
		float a2 = a_2;
		float b1 = b_1;
		float b2 = b_2;
		float c0 = 1.0;
		float d0 = 0.0;

		std::array<float, 7> coeffs{ a0, a1, a2, b1, b2, c0, d0 };
		m_biquad.setCoeffs(coeffs);

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