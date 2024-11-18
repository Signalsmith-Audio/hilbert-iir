/* Signalsmith's Hilbert IIR: A single-file, dependency-free Hilbert filter.

Copyright (c) 2024 Geraint Luff / Signalsmith Audio Ltd.

Released under 0BSD: BSD Zero Clause License
*/
#ifndef SIGNALSMITH_HILBERT_IIR_H
#define SIGNALSMITH_HILBERT_IIR_H

#include <complex>
#include <array>
#include <vector>

namespace signalsmith { namespace hilbert {

// produced with design.py
template<typename Sample>
struct HilbertIIRCoeffs {
	static constexpr int order = 12;
	std::array<std::complex<Sample>, order> coeffs{{
		{Sample(-0.000224352093802), Sample(0.00543499018201)},
		{Sample(0.0107500557815), Sample(-0.0173890685681)},
		{Sample(-0.0456795873917), Sample(0.0229166931429)},
		{Sample(0.11282500582), Sample(0.00278413661237)},
		{Sample(-0.208067578452), Sample(-0.104628958675)},
		{Sample(0.28717837501), Sample(0.33619239719)},
		{Sample(-0.254675294431), Sample(-0.683033899655)},
		{Sample(0.0481081835026), Sample(0.954061589374)},
		{Sample(0.227861357867), Sample(-0.891273574569)},
		{Sample(-0.365411839137), Sample(0.525088317271)},
		{Sample(0.280729061131), Sample(-0.155131206606)},
		{Sample(-0.0935061787728), Sample(0.00512245855404)}
	}};
	std::array<std::complex<Sample>, order> poles{{
		{Sample(-0.00495335976478), Sample(0.0092579876872)},
		{Sample(-0.017859491302), Sample(0.0273493725543)},
		{Sample(-0.0413714373155), Sample(0.0744756910287)},
		{Sample(-0.0882148408885), Sample(0.178349677457)},
		{Sample(-0.17922965812), Sample(0.39601340223)},
		{Sample(-0.338261800753), Sample(0.829229533354)},
		{Sample(-0.557688699732), Sample(1.61298538328)},
		{Sample(-0.735157736148), Sample(2.79987398682)},
		{Sample(-0.719057381172), Sample(4.16396166128)},
		{Sample(-0.517871025209), Sample(5.29724826804)},
		{Sample(-0.280197469471), Sample(5.99598602388)},
		{Sample(-0.0852751354531), Sample(6.3048492377)}
	}};
	Sample direct = 0.000262057212648;
};

template<typename Sample>
struct HilbertIIR {
	using Complex = std::complex<Sample>;
	static constexpr int order = HilbertIIRCoeffs<Sample>::order;
	
	HilbertIIR(Sample sampleRate=48000, int channels=1, Sample passbandGain=2) {
		HilbertIIRCoeffs<Sample> coeffs;
		
		Sample freqFactor = std::min<Sample>(0.46, 20000/sampleRate);
		direct = coeffs.direct*2*passbandGain*freqFactor;
		for (int i = 0; i < order; ++i) {
			Complex coeff = coeffs.coeffs[i]*freqFactor*passbandGain;
			coeffsR[i] = coeff.real();
			coeffsI[i] = coeff.imag();
			Complex pole = std::exp(coeffs.poles[i]*freqFactor);
			polesR[i] = pole.real();
			polesI[i] = pole.imag();
		}
		states.resize(channels);
		reset();
	}
	
	void reset() {
		for (auto &state : states) {
			for (auto &v : state.real) v = 0;
			for (auto &v : state.imag) v = 0;
		}
	}
	
	Complex operator()(Sample x, int channel=0) {
		// Really we're just doing: state[i] = state[i]*poles[i] + x*coeffs[i]
		// but std::complex is slow without -ffast-math, so we've unwrapped it
		State state = states[channel], newState;
		for (int i = 0; i < order; ++i) {
			newState.real[i] = state.real[i]*polesR[i] - state.imag[i]*polesI[i] + x*coeffsR[i];
		}
		for (int i = 0; i < order; ++i) {
			newState.imag[i] = state.real[i]*polesI[i] + state.imag[i]*polesR[i] + x*coeffsI[i];
		}
		states[channel] = newState;

		Sample resultR = x*direct;
		for (int i = 0; i < order; ++i) {
			resultR += newState.real[i];
		}
		Sample resultI = 0;
		for (int i = 0; i < order; ++i) {
			resultI += newState.imag[i];
		}
		return {resultR, resultI};
	}

	Complex operator()(Complex x, int channel=0) {
		State state = states[channel], newState;
		for (int i = 0; i < order; ++i) {
			newState.real[i] = state.real[i]*polesR[i] - state.imag[i]*polesI[i] + x.real()*coeffsR[i] - x.imag()*coeffsI[i];
		}
		for (int i = 0; i < order; ++i) {
			newState.imag[i] = state.real[i]*polesI[i] + state.imag[i]*polesR[i] + x.real()*coeffsI[i] + x.imag()*coeffsR[i];
		}
		states[channel] = newState;

		Sample resultR = x.real()*direct;
		for (int i = 0; i < order; ++i) {
			resultR += newState.real[i];
		}
		Sample resultI = x.imag()*direct;;
		for (int i = 0; i < order; ++i) {
			resultI += newState.imag[i];
		}
		return {resultR, resultI};
	}
private:
	using Array = std::array<Sample, order>;
	Array coeffsR, coeffsI, polesR, polesI;
	struct State {
		Array real, imag;
	};
	std::vector<State> states;
	Sample direct;
};

}} // namespace
#endif // include guard
