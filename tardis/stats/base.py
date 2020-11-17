import numpy as np


def get_trivial_poisson_uncertainty(model):
    """"""
    emitted_nu = model.montecarlo_nu[model.montecarlo_luminosity >= 0]
    emitted_luminosity = model.montecarlo_luminosity[
        model.montecarlo_luminosity >= 0
    ]
    freq_bins = model.tardis_config.spectrum.frequency.value

    bin_counts = np.histogram(emitted_nu, bins=freq_bins)[0]

    uncertainty = np.sqrt(bin_counts) * np.mean(emitted_luminosity)

    return uncertainty / (freq_bins[1] - freq_bins[0])
