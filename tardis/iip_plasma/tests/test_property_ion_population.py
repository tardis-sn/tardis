import numpy as np
import pytest

pytestmark = pytest.mark.skip("Skipping tests due to old format")


def test_phi_saha_lte(beta_rad, g_electron, ionization_data, phi_saha_lte):
    assert phi_saha_lte.shape == (2, 20)
    assert np.all(
        np.isclose(
            g_electron * 4 * np.exp(-ionization_data.loc[2, 1] * beta_rad),
            phi_saha_lte.loc[2].loc[1],
        )
    )


def test_ion_number_density(
    phi_saha_lte, ion_number_density, electron_densities
):
    ion_numbers = ion_number_density.index.get_level_values(1).values
    ion_numbers = ion_numbers.reshape((ion_numbers.shape[0], 1))
    n_electron_from_ions = (ion_number_density.values * ion_numbers).sum(axis=0)
    assert (
        np.all(
            np.abs(
                (n_electron_from_ions - electron_densities) / electron_densities
            )
            < 0.05
        )
        == True
    )


def test_electron_densities(electron_densities):
    assert np.allclose(electron_densities, 1.181197e09)


def test_radiation_field_correction(delta):
    assert np.allclose(delta.loc[2].loc[2], 0.000807200897)


def test_phi_saha_nebular(phi_saha_nebular):
    assert np.allclose(phi_saha_nebular.loc[2].loc[1], 4.496818e08)
