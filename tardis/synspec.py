#creation of spectra at the end of the run
import numpy as np
import constants
#deprectated will use astropy in the end
import pyspec

def get_lambda_spec(nu, energy, wave_min, wave_max, samples=50):
    wave = 1e8 * constants.c / np.array(nu)
    nu_bins = np.linspace(constants.c / wave_min, constants.c / wave_max, samples)
    wavelength_bins = 1e8 * constants.c / np.array(nu_bins)
    wavelength = wavelength_bins[:-1] + np.diff(wavelength_bins)*0.5
    intens = np.histogram(wave, weights=energy, bins=wavelength_bins)[0]
    
    return wavelength, intens
    