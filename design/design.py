"""
Designs a continous-time elliptic filter, shifts the results to make a Hilbert, and writing the results as a C++ struct
"""

import numpy
import scipy
import scipy.signal
import scipy.optimize

filter_order = 12
filter_ripple_db = 0.5
filter_stop_db = 90
filter_squinch = 0.5
filter_offset = 0.01; # moves the lower cutoff slightly above 0

#------- Create the filter

(filter_z, filter_p, filter_k) = scipy.signal.ellip(filter_order, filter_ripple_db, filter_stop_db, numpy.pi, analog=True, output='zpk')

# Squinch (improves DC transition at the expense of upper band)
filter_z -= numpy.pi*1j
filter_p -= numpy.pi*1j
filter_z /= 1 + filter_z*1j*filter_squinch
filter_p /= 1 + filter_p*1j*filter_squinch
K = (1 + numpy.pi*2*filter_squinch)
filter_z = filter_z*K + numpy.pi*2j + filter_offset*1j
filter_p = filter_p*K + numpy.pi*2j + filter_offset*1j

(filter_b, filter_a) = scipy.signal.zpk2tf(filter_z, filter_p, filter_k)

residue_coeffs, residue_poles, residue_direct = scipy.signal.residue(filter_b, filter_a)
residue_coeffs *= 0.5
residue_direct *= 2*K

# Write the results out as a C++ class
print("""template<typename Sample>
struct HilbertIIRCoeffs {
	static constexpr int order = %i;
	std::array<std::complex<Sample>, order> coeffs{{
		%s
	}};
	std::array<std::complex<Sample>, order> poles{{
		%s
	}};
	Sample direct = %s;
};"""%(
	filter_order,
	",\n\t\t".join(["{Sample(%s), Sample(%s)}"%((2*p.real).astype(str), (2*p.imag).astype(str)) for p in residue_coeffs]),
	",\n\t\t".join(["{Sample(%s), Sample(%s)}"%(p.real.astype(str), p.imag.astype(str)) for p in residue_poles]),
	residue_direct[0].real
))
