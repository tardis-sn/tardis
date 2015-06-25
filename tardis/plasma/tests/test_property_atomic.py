import numpy as np

def test_levels_property(levels):
    assert np.isclose(levels.ix[2].ix[0].ix[1]['energy'], 3.17545416e-11)

def test_lines_property(lines):
    assert np.isclose(lines.ix[564954]['wavelength'], 10833.307)
    assert lines.index[124] == 564954

def test_lines_lower_level_index_property(lines_lower_level_index):
    assert lines_lower_level_index[9] == 0

def test_lines_upper_level_index_property(lines_upper_level_index):
    assert lines_upper_level_index[9] == 30

def test_atomic_mass_property(atomic_mass):
    assert np.isclose(atomic_mass.ix[2], 6.6464755216973998e-24)

def test_ionization_data_property(ionization_data):
    assert np.isclose(float(ionization_data.ix[2].ix[1]), 3.9393336e-11)

def test_zeta_data_property(zeta_data):
    assert zeta_data.shape == (3, 20)
    assert np.isclose(zeta_data.ix[2].ix[1][2000], 0.4012)