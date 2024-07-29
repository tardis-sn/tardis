import astropy.units as u
import numpy as np


def calculate_emitted_luminosity(
    emitted_packet_nu,
    emitted_packet_luminosity,
    luminosity_nu_start=0 * u.Hz,
    luminosity_nu_end=np.inf * u.Hz,
):
    """
    Calculate emitted luminosity.

    Parameters
    ----------
    emitted_packet_nu :
    emitted_packet_luminosity :
    luminosity_nu_start : astropy.units.Quantity
    luminosity_nu_end : astropy.units.Quantity

    Returns
    -------
    astropy.units.Quantity
    """
    luminosity_wavelength_filter = (emitted_packet_nu > luminosity_nu_start) & (
        emitted_packet_nu < luminosity_nu_end
    )

    return emitted_packet_luminosity[luminosity_wavelength_filter].sum()


def calculate_reabsorbed_luminosity(
    reabsorbed_packet_nu,
    reabsorbed_packet_luminosity,
    luminosity_nu_start=0 * u.Hz,
    luminosity_nu_end=np.inf * u.Hz,
):
    """
    Calculate reabsorbed luminosity.

    Parameters
    ----------
    reabsorbed_packet_nu :
    reabsorbed_packet_luminosity :
    luminosity_nu_start : astropy.units.Quantity
    luminosity_nu_end : astropy.units.Quantity

    Returns
    -------
    astropy.units.Quantity
    """
    luminosity_wavelength_filter = (
        reabsorbed_packet_nu > luminosity_nu_start
    ) & (reabsorbed_packet_nu < luminosity_nu_end)

    return reabsorbed_packet_luminosity[luminosity_wavelength_filter].sum()
