import pytest
import numpy.testing as npt

import tardis.energy_input.calculate_opacity as calculate_opacity


@pytest.mark.parametrize(
    ["electron_number_density", "energy", "expected"],
    [
        (1.0, 511.0, 0.0856559215396578),
        (1e-2, 255.5, 0.0011191738319133954),
        (1e5, 511.0e7, 0.012909632812042571),
    ],
)
def test_compton_opacity_calculation(energy, electron_number_density, expected):
    """
    Parameters
    ----------
    energy : float
    electron_number_density : float
    """
    opacity = calculate_opacity.compton_opacity_calculation(
        energy, electron_number_density
    )

    npt.assert_almost_equal(opacity, expected)


@pytest.mark.parametrize(
    ["ejecta_density", "energy", "iron_group_fraction", "expected"],
    [
        (1.0, 511.0, 0.0, 0.0001497022737728184),
        (1e-2, 255.5, 0.5, 8.903267700390038e-05),
        (1e-2, 255.5, 0.25, 5.1069068712110425e-05),
        (1e5, 511.0e7, 1.0, 0.0),
    ],
)
def test_photoabsorption_opacity_calculation(
    energy, ejecta_density, iron_group_fraction, expected
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

    npt.assert_almost_equal(opacity, expected)


@pytest.mark.parametrize(
    ["ejecta_density", "energy", "iron_group_fraction", "expected"],
    [
        (1.0, 511.0, 0.0, 0.0),
        (1e-2, 1500, 0.5, 2.743980356831218e-06),
        (1e-2, 1200, 0.25, 8.846018943383742e-06),
        (1e5, 511.0e7, 1.0, 1113145501.6992927),
    ],
)
def test_pair_creation_opacity_calculation(
    energy, ejecta_density, iron_group_fraction, expected
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

    npt.assert_almost_equal(opacity, expected)
