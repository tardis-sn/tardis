import numpy as np
import pytest

pytestmark = pytest.mark.skip("Skipping tests due to old format")


def test_levels_property(excitation_energy):
    assert np.isclose(excitation_energy.loc[2].loc[0].loc[1], 3.17545416e-11)


def test_lines_property(lines):
    assert np.isclose(lines.loc[564954]["wavelength"], 10833.307)
    assert lines.index[124] == 564954


def test_lines_lower_level_index_property(lines_lower_level_index):
    assert lines_lower_level_index[9] == 0


def test_lines_upper_level_index_property(lines_upper_level_index):
    assert lines_upper_level_index[9] == 30


def test_atomic_mass_property(atomic_mass):
    assert np.isclose(atomic_mass.loc[2], 6.6464755216973998e-24)


def test_ionization_data_property(ionization_data):
    assert np.isclose(float(ionization_data.loc[2].loc[1]), 3.9393336e-11)


def test_zeta_data_property(zeta_data):
    assert np.isclose(zeta_data.loc[2].loc[1][2000], 0.4012)
