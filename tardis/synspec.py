#creation of spectra at the end of the run
import numpy as np
import constants
#deprectated will use astropy in the end
from pyspec import oned

def get_lambda_spec(nu, energy, wave_min, wave_max, samples=50):
    nu = nu[nu >= 0]
    energy = energy[energy >= 0]
    wave = 1e8 * constants.c / np.array(nu)
    nu_bins = np.linspace(constants.c / wave_min, constants.c / wave_max, samples)
    wavelength_bins = 1e8 * constants.c / (np.array(nu_bins))
    wavelength = wavelength_bins[:-1] + np.diff(wavelength_bins) * 0.5
    intens = np.histogram(wave, weights=energy, bins=wavelength_bins)[0]

    return oned.onedspec(wavelength, intens * constants.c / wavelength ** 2, mode='waveflux')


class tardis_result(object):
    def __init__(self,
                 nu_input,
                 energy_of_packet,
                 nu,
                 energy,
                 nu_reabsorbed,
                 energy_reabsorbed,
                 track_t_rads,
                 track_ws,
                 track_t_inner,
                 wave_min=500,
                 wave_max=20000,
                 samples=2000):
        wave, intens = get_lambda_spec(nu_input, np.ones_like(nu_input) * energy_of_packet, wave_min * 1e-8,
            wave_max * 1e-8, samples=samples)
        self.input_spec = oned.onedspec(wave, intens, mode='waveflux')

        nu_filter_emit = nu != 0.0
        nu = nu[nu_filter_emit]
        energy = energy[nu_filter_emit]
        wave, intens = get_lambda_spec(nu, energy, wave_min * 1e-8, wave_max * 1e-8, samples=samples)
        self.output_spec = oned.onedspec(wave, intens, mode='waveflux')

        nu_filter = nu_reabsorbed != 0.0
        nu_reabsorbed = nu_reabsorbed[nu_filter]
        energy_reabsorbed = energy_reabsorbed[nu_filter]

        wave, intens = get_lambda_spec(nu_reabsorbed, energy_reabsorbed, wave_min * 1e-8, wave_max * 1e-8,
            samples=samples)
        self.reabsorbed_spec = oned.onedspec(wave, intens, mode='waveflux')

