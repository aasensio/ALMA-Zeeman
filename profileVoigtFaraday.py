from scipy.special import wofz
#-----------------------------------------
# Returns the Voigt function for an axis of wavelengths l and damping parameter a
#-----------------------------------------
def profileVoigtFaraday(l, a):
	z = l + 1j*a
	return wofz(z).real, wofz(z).imag