import pytest
import numpy.testing as npt

import tardis.energy_input.calculate_opacity as calculate_opacity


@pytest.mark.parametrize(
    ["electron_number_density", "energy", "expected"],
    [
        (1.0e11, 511.0, 2.865396624016367e-14),
        (1e15, 255.5, 3.743906253489761e-10),
        (1e5, 511.0e7, 4.318577913631238e-26),
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
        (1.0, 511.0, 0.0, 0.00015028056615643418),
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
        (1e5, 511.0e7, 1.0, 1111355719.7411418),
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
