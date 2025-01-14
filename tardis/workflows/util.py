import numpy as np
from astropy import constants as const


import numpy as np
from astropy import constants as const, units as u


def get_tau_integ(plasma, opacity_state, simulation_state, bin_size=10):
    """Estimate the integrated mean optical depth at each velocity bin

    Parameters
    ----------
    plasma : tardis.plasma.BasePlasma
        The tardis legacy plasma
    simulation_state : tardis.model.base.SimulationState
        the current simulation state
    bin_size : int, optional.  Default : 10
        bin size for the aggregation of line opacities

    Returns
    -------
    dict
        rosseland : np.ndarray
            Roassland Mean Optical Depth
        planck : np.ndarray
            Planck Mean Optical Depth
    """
    index = plasma.atomic_data.lines.nu.index
    freqs = plasma.atomic_data.lines.nu.values * u.Hz
    order = np.argsort(freqs)
    freqs = freqs[order]
    taus = opacity_state.tau_sobolev.values[order]

    check_bin_size = True
    while check_bin_size:
        extra = bin_size - len(freqs) % bin_size
        extra_freqs = (np.arange(extra + 1) + 1) * u.Hz
        extra_taus = np.zeros((extra + 1, taus.shape[1]))
        freqs = np.hstack((extra_freqs, freqs))
        taus = np.vstack((extra_taus, taus))

        bins_low = freqs[:-bin_size:bin_size]
        bins_high = freqs[bin_size::bin_size]
        delta_nu = bins_high - bins_low
        n_bins = len(delta_nu)

        if np.any(delta_nu == 0):
            bin_size += 2
        else:
            check_bin_size = False

    taus = taus[1 : n_bins * bin_size + 1]
    freqs = freqs[1 : n_bins * bin_size + 1]

    ct = simulation_state.time_explosion * const.c
    t_rad = simulation_state.radiation_field_state.temperature

    def B(nu, T):
        return (
            2
            * const.h
            * nu**3
            / const.c**2
            / (np.exp(const.h * nu / (const.k_B * T)) - 1)
        )

    def U(nu, T):
        return (
            B(nu, T) ** 2 * (const.c / nu) ** 2 * (2 * const.k_B * T**2) ** -1
        )

    kappa_exp = (
        (bins_low / delta_nu).reshape(-1, 1)
        / ct
        * (1 - np.exp(-taus.reshape(n_bins, bin_size, -1))).sum(axis=1)
    )
    kappa_thom = plasma.electron_densities.values * u.cm ** (-3) * const.sigma_T
    Bdnu = B(bins_low.reshape(-1, 1), t_rad.reshape(1, -1)) * delta_nu.reshape(
        -1, 1
    )
    kappa_planck = kappa_thom + (Bdnu * kappa_exp).sum(axis=0) / (
        Bdnu.sum(axis=0)
    )

    udnu = U(bins_low.reshape(-1, 1), t_rad.reshape(1, -1)) * delta_nu.reshape(
        -1, 1
    )
    kappa_tot = kappa_thom + kappa_exp
    kappa_rosseland = (
        (udnu * kappa_tot**-1).sum(axis=0) / (udnu.sum(axis=0))
    ) ** -1

    dr = simulation_state.geometry.r_outer - simulation_state.geometry.r_inner
    dtau = kappa_planck * dr
    planck_integ_tau = np.cumsum(dtau[::-1])[::-1]
    rosseland_integ_tau = np.cumsum((kappa_rosseland * dr)[::-1])[::-1]

    return {"rosseland": rosseland_integ_tau, "planck": planck_integ_tau}
