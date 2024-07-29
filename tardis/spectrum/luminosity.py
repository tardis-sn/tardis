import astropy.units as u
import numpy as np


def calculate_filtered_luminosity(
    packet_nu,
    packet_luminosity,
    luminosity_nu_start=0 * u.Hz,
    luminosity_nu_end=np.inf * u.Hz,
):
    """
    Calculate total luminosity within a filter range.

    Parameters
    ----------
    packet_nu : astropy.units.Quantity
    packet_luminosity : astropy.units.Quantity
    luminosity_nu_start : astropy.units.Quantity
    luminosity_nu_end : astropy.units.Quantity

    Returns
    -------
    astropy.units.Quantity
    """
    luminosity_wavelength_filter = (packet_nu > luminosity_nu_start) & (
        packet_nu < luminosity_nu_end
    )

    return packet_luminosity[luminosity_wavelength_filter].sum()
