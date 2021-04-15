import pytest
import astropy.units as u
import numpy.testing as npt
import numpy as np

import tardis.energy_input.calculate_opacity as calculate_opacity
import tardis.energy_input.util as util
from tardis import constants as const

# this is terrible practice
MASS_SI = 28.085 * const.m_p.to(u.g).value
MASS_FE = 55.845 * const.m_p.to(u.g).value
M_P = const.m_p.to(u.g).value


@pytest.mark.parametrize(
    ["electron_number_density", "energy"],
    [(1.0, 511.0), (1e-2, 255.5), (1e5, 511.0e7)],
)
def test_compton_opacity_calculation(energy, electron_number_density):
    """
    Parameters
    ----------
    energy : float
    electron_number_density : float
    """
    opacity = calculate_opacity.compton_opacity_calculation(
        energy, electron_number_density
    )

    kappa = util.kappa_calculation(energy)

    sigma_T = const.sigma_T.cgs.value

    a = 1.0 + 2.0 * kappa

    sigma_KN = (
        3.0
        / 4.0
        * sigma_T
        * (
            (1.0 + kappa)
            / kappa ** 3.0
            * ((2.0 * kappa * (1.0 + kappa)) / a - np.log(a))
            + 1.0 / (2.0 * kappa) * np.log(a)
            - (1.0 + 3 * kappa) / a ** 2.0
        )
    )

    expected = electron_number_density / (M_P * 2) * sigma_KN

    npt.assert_almost_equal(opacity, expected)


@pytest.mark.parametrize(
    ["ejecta_density", "energy", "iron_group_fraction"],
    [
        (1.0, 511.0, 0.0),
        (1e-2, 255.5, 0.5),
        (1e-2, 255.5, 0.25),
        (1e5, 511.0e7, 1.0),
    ],
)
def test_photoabsorption_opacity_calculation(
    energy, ejecta_density, iron_group_fraction
):
    """
    Parameters
    ----------
    energy : float
    ejecta_density : float
    iron_group_fraction : float
    """
    opacity = calculate_opacity.photoabsorption_opacity_calculation(
        energy, ejecta_density, iron_group_fraction
    )

    Si_opacity = (
        1.16e-24
        * (energy / 100.0) ** -3.13
        * ejecta_density
        / MASS_SI
        * (1.0 - iron_group_fraction)
    )

    Fe_opacity = (
        25.7e-24
        * (energy / 100.0) ** -3.0
        * ejecta_density
        / MASS_FE
        * (1.0 - iron_group_fraction)
    )

    expected = Si_opacity + Fe_opacity

    npt.assert_almost_equal(opacity, expected)


@pytest.mark.parametrize(
    ["ejecta_density", "energy", "iron_group_fraction"],
    [
        (1.0, 511.0, 0.0),
        (1e-2, 255.5, 0.5),
        (1e-2, 255.5, 0.25),
        (1e5, 511.0e7, 1.0),
    ],
)
def test_pair_creation_opacity_calculation(
    energy, ejecta_density, iron_group_fraction
):
    """
    Parameters
    ----------
    energy : float
    ejecta_density : float
    iron_group_fraction : float
    """
    opacity = calculate_opacity.pair_creation_opacity_calculation(
        energy, ejecta_density, iron_group_fraction
    )

    Z_Si = 14
    Z_Fe = 26

    Si_proton_ratio = Z_Si ** 2.0 / MASS_SI
    Fe_proton_ratio = Z_Fe ** 2.0 / MASS_FE

    multiplier = ejecta_density * (
        Si_proton_ratio * (1.0 - iron_group_fraction)
        + Fe_proton_ratio * iron_group_fraction
    )

    if energy > 1.022e3 and energy < 1.5e3:
        expected = multiplier * 1.0063 * (energy / 1.0e3 - 1.022) * 1.0e-27
    else:
        expected = (
            multiplier * (0.0481 + 0.301 * (energy / 1.0e3 - 1.5)) * 1.0e-27
        )

    if expected < 0.0:
        expected = 0.0

    npt.assert_almost_equal(opacity, expected)
