# Hilbert IIR

This is a (short, dependency-free, C++ header-only) IIR Hilbert filter with an [0BSD License](LICENSE.txt).

It uses a pre-baked continuous-time filter, which is converted into a discrete-time filter (with an upper cutoff of 20kHz) for a given sample-rate.  This gives it the same phase/amplitude responses at higher sample-rates (as per my API London 2024 talk).

## How to use it:

```cpp
#import "hilbert.h"
using HilbertIIR = signalsmith::hilbert::HilbertIIR<float>;
```

Single-channel:

```cpp
HilbertIIR hilbert(48000); // 48kHz

float in = /*...*/;
std::complex<float> out = hilbert(in);
```

Multi-channel:

```cpp
HilbertIIR hilbert(96000, 2); // 96kHz stereo

float in0 = /*...*/, in1 = /*...*/;
std::complex<float> out0 = hilbert(in0, 0);
std::complex<float> out1 = hilbert(in1, 1);
```

Reset everything to 0:

```cpp
hilbert.reset();
```

## Response

Here's the impulse/spectrum when running at 44.1kHz:

<img src="design/plots/plot-44100@2x.png" width="472" height="422">

Here's the impulse/spectrum when running at 96kHz:

<img src="design/plots/plot-96000@2x.png" width="472" height="422">

### Filter design / IIT

A continuous-time filter is designed in Python / SciPy (`design/design.py`).  It's [rephrased as](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.residue.html) the sum of parallel complex 1-pole filters, and printed as C++ code.

The C++ code (`hilbert.h`) then converts this to a discrete-time filter using the [Impulse Invariant Transform](https://en.wikipedia.org/wiki/Impulse_invariance#Effect_on_poles_in_system_function), so we now have the sum of parallel *discrete-time* complex 1-pole filters.  This is what we use to compute the result.

### Plots

The code in `design/main.cpp` generates SVG plots of the impulse/spectrum, so if you change the filter design you can see the results.

## License

The license is currently [0BSD](LICENSE.txt), but if you need anything else, get in touch.
