#ifndef SIGNALSMITH_IIR_HILBERT_H
#define SIGNALSMITH_IIR_HILBERT_H

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
		{Sample(0.0481081835038), Sample(0.954061589374)},
		{Sample(0.227861357868), Sample(-0.891273574562)},
		{Sample(-0.365411839149), Sample(0.525088317279)},
		{Sample(0.280729061125), Sample(-0.155131206624)},
		{Sample(-0.0935061787568), Sample(0.00512245855511)}
	}};
	std::array<std::complex<Sample>, order> poles{{
		{Sample(-0.00495335976478), Sample(-0.000742012312798)},
		{Sample(-0.017859491302), Sample(0.0173493725543)},
		{Sample(-0.0413714373155), Sample(0.0644756910287)},
		{Sample(-0.0882148408885), Sample(0.168349677457)},
		{Sample(-0.17922965812), Sample(0.38601340223)},
		{Sample(-0.338261800753), Sample(0.819229533354)},
		{Sample(-0.557688699732), Sample(1.60298538328)},
		{Sample(-0.735157736147), Sample(2.78987398682)},
		{Sample(-0.719057381177), Sample(4.15396166127)},
		{Sample(-0.517871025196), Sample(5.28724826806)},
		{Sample(-0.280197469484), Sample(5.98598602386)},
		{Sample(-0.0852751354486), Sample(6.29484923772)}
	}};
	Sample direct = 3.16372510007e-05;
};

template<typename Sample>
struct HilbertIIR {
	using Complex = std::complex<Sample>;
	static constexpr int order = HilbertIIRCoeffs<Sample>::order;
	
	HilbertIIR(Sample sampleRate=48000, int channels=1) {
		Sample freqFactor = std::min<Sample>(0.5, 20000/sampleRate);
		for (auto &z : coeffs.coeffs) {
			z *= freqFactor;
		}
		for (auto &z : coeffs.poles) {
			z = std::exp(freqFactor*z);
		}
		states.resize(channels);
		reset();
	}
	
	void reset() {
		for (auto &state : states) {
			for (auto &z : state) z = 0;
		}
	}
	
	Complex operator()(Sample x, int channel=0) {
		auto &state = states[channel];
		for (int i = 0; i < order; ++i) {
			state[i] = state[i]*coeffs.poles[i] + x*coeffs.coeffs[i];
		}
		Complex result = x*coeffs.direct;
		for (int i = 0; i < order; ++i) {
			result += state[i];
		}
		return result;
	}
private:
	HilbertIIRCoeffs<Sample> coeffs;
	using State = std::array<Complex, order>;
	std::vector<State> states;
};

}} // namespace
#endif // include guard
