#include "../hilbert.h"

#include "./plot.h"
#include "./simple-fft.h"

using Sample = float;
using Complex = std::complex<Sample>;

void plot(Sample sampleRate=48000) {
	int impulseLength = 100, fftSize = 32768;
	signalsmith::hilbert::HilbertIIR<Sample> hilbert(sampleRate, 1, 1);
	
	signalsmith::plot::Figure figure;
	auto &impulsePlot = figure(0, 0).plot(400, 80);
	impulsePlot.x.major(0).minors(25, 50, 75, 100);
	impulsePlot.y.major(0);
	auto &impulseReal = impulsePlot.line();
	auto &impulseImag = impulsePlot.line();
	
	std::vector<Complex> impulse(fftSize), spectrum(fftSize);
	
	for (int i = 0; i < fftSize; ++i) {
		Sample x = (i == 0);
		Complex y = hilbert(x);
		if (i < impulseLength) {
			impulseReal.add(i, y.real());
			impulseImag.add(i, y.imag());
		}
		
		impulse[i] = y;
	}
	
	signalsmith::fft2::SimpleFFT<Sample> fft(fftSize);
	fft.fft(fftSize, impulse.data(), spectrum.data());

	auto &spectrumPlot = figure(0, 1).plot(400, 100);
	auto &spectrumPlotZoom = figure(0, 2).plot(400, 100);
	spectrumPlot.x.major(0).major(sampleRate/2, "N").major(-sampleRate/2, "-N").minor(20000, "20k").label("Hz");
	spectrumPlot.y.major(0).minors(-120, -90, -60, -30).label("dB");
	spectrumPlotZoom.x.major(0).minors(-200, 200).linear(-400, 400).label("Hz");
	spectrumPlotZoom.y.major(0).minors(-60, -40, -20).linear(-80, 3).label("dB");
	auto &dbLine = spectrumPlot.line();
	auto &dbLine2 = spectrumPlotZoom.line();
	for (int i = 0; i < fftSize; ++i) {
		int b = (i + fftSize/2)%fftSize;
		Sample f = (i - fftSize/2)*sampleRate/fftSize;
		Complex bin = spectrum[b];
		Sample db = std::max<Sample>(-120, 10*std::log10(std::norm(bin) + 1e-30));
		dbLine.add(f, db);
		dbLine2.add(f, db);
	}

	figure.write("plot-" + std::to_string(int(sampleRate)) + ".svg");
}
int main() {
	plot(44100);
	plot(48000);
	plot(96000);
	plot(192000);
}
