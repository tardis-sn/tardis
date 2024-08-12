import astropy.units as u
import numpy as np
import pytest

from tardis.spectrum.luminosity import (
    calculate_filtered_luminosity,
)


@pytest.mark.parametrize(
    "packet_nu, packet_luminosity, luminosity_nu_start, luminosity_nu_end, expected",
    [
        # All frequencies within the range
        (
            np.array([1, 2, 3]) * u.Hz,
            np.array([10, 20, 30]) * u.erg / u.s,
            0 * u.Hz,
            4 * u.Hz,
            60 * u.erg / u.s,
        ),
        # All frequencies outside the range
        (
            np.array([1, 2, 3]) * u.Hz,
            np.array([10, 20, 30]) * u.erg / u.s,
            4 * u.Hz,
            5 * u.Hz,
            0 * u.erg / u.s,
        ),
        # Mix of frequencies within and outside the range
        (
            np.array([1, 2, 3, 4]) * u.Hz,
            np.array([10, 20, 30, 40]) * u.erg / u.s,
            2 * u.Hz,
            4 * u.Hz,
            30 * u.erg / u.s,
        ),
        # Edge case: Frequencies exactly on the boundary
        (
            np.array([1, 2, 3, 4]) * u.Hz,
            np.array([10, 20, 30, 40]) * u.erg / u.s,
            2 * u.Hz,
            3 * u.Hz,
            0 * u.erg / u.s,
        ),
    ],
)
def test_calculate_filtered_luminosity(
    packet_nu,
    packet_luminosity,
    luminosity_nu_start,
    luminosity_nu_end,
    expected,
):
    result = calculate_filtered_luminosity(
        packet_nu,
        packet_luminosity,
        luminosity_nu_start,
        luminosity_nu_end,
    )
    assert result == expected
