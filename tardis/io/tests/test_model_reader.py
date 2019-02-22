import os

from astropy import units as u
import numpy as np
import pytest

import tardis
from tardis.io.config_reader import Configuration
from tardis.io.model_reader import (
    read_artis_density, read_simple_ascii_abundances, read_csv_composition, read_uniform_abundances, read_cmfgen_density, read_cmfgen_composition)

data_path = os.path.join(tardis.__path__[0], 'io', 'tests', 'data')

@pytest.fixture
def artis_density_fname():
    return os.path.join(data_path, 'artis_model.dat')

@pytest.fixture
def artis_abundances_fname():
    return os.path.join(data_path, 'artis_abundances.dat')

@pytest.fixture
def cmfgen_fname():
    return os.path.join(data_path, 'cmfgen_model.csv')

@pytest.fixture
def csv_composition_fname():
    return os.path.join(data_path, 'csv_composition.csv')


@pytest.fixture
def isotope_uniform_abundance():
    config_path = os.path.join(
        data_path, 'tardis_configv1_isotope_uniabund.yml')
    config = Configuration.from_yaml(config_path)
    return config.model.abundances

def test_simple_read_artis_density(artis_density_fname):
    time_of_model, velocity, mean_density = read_artis_density(artis_density_fname)

    assert np.isclose(0.00114661 * u.day, time_of_model, atol=1e-7 * u.day)
    assert np.isclose(mean_density[23], 0.2250048 * u.g / u.cm**3, atol=1.e-6
        * u.g / u.cm**3)
    assert len(mean_density) == 69
    assert len(velocity) == len(mean_density) + 1

# Artis files are currently read with read ascii files function
def test_read_simple_ascii_abundances(artis_abundances_fname):
    index, abundances = read_simple_ascii_abundances(artis_abundances_fname)
    assert len(abundances.columns) == 69
    assert np.isclose(abundances[23].loc[2], 2.672351e-08 , atol=1.e-12)


def test_read_simple_isotope_abundances(csv_composition_fname):
    index, abundances, isotope_abundance = read_csv_composition(
        csv_composition_fname)
    assert np.isclose(abundances.loc[6, 8], 0.5, atol=1.e-12)
    assert np.isclose(abundances.loc[12, 5], 0.8, atol=1.e-12)
    assert np.isclose(abundances.loc[14, 1], 0.1, atol=1.e-12)
    assert np.isclose(isotope_abundance.loc[(28, 56), 0], 0.4, atol=1.e-12)
    assert np.isclose(isotope_abundance.loc[(28, 58), 2], 0.7, atol=1.e-12)
    assert abundances.shape == (4, 10)
    assert isotope_abundance.shape == (2, 10)


def test_read_cmfgen_isotope_abundances(cmfgen_fname):
    index, abundances, isotope_abundance = read_cmfgen_composition(
        cmfgen_fname)
    assert np.isclose(abundances.loc[6, 8], 0.5, atol=1.e-12)
    assert np.isclose(abundances.loc[12, 5], 0.8, atol=1.e-12)
    assert np.isclose(abundances.loc[14, 1], 0.3, atol=1.e-12)
    assert np.isclose(isotope_abundance.loc[(28, 56), 0], 0.5, atol=1.e-12)
    assert np.isclose(isotope_abundance.loc[(28, 58), 1], 0.7, atol=1.e-12)
    assert abundances.shape == (4, 9)
    assert isotope_abundance.shape == (2, 9)


def test_read_uniform_abundances(isotope_uniform_abundance):
    abundances, isotope_abundance = read_uniform_abundances(
        isotope_uniform_abundance, 20)
    assert np.isclose(abundances.loc[8, 2], 0.19, atol=1.e-12)
    assert np.isclose(abundances.loc[20, 5], 0.03, atol=1.e-12)
    assert np.isclose(isotope_abundance.loc[(28, 56), 15], 0.05, atol=1.e-12)
    assert np.isclose(isotope_abundance.loc[(28, 58), 2], 0.05, atol=1.e-12)


def test_simple_read_cmfgen_density(cmfgen_fname):
    time_of_model, velocity, mean_density, electron_densities, temperature = read_cmfgen_density(
        cmfgen_fname)

    assert np.isclose(0.976 * u.day, time_of_model, atol=1e-7 * u.day)
    assert np.isclose(mean_density[4], 4.2539537e-09 * u.g / u.cm**3, atol=1.e-6
                      * u.g / u.cm**3)
    assert np.isclose(electron_densities[5], 2.6e+14 * u.cm**-3, atol=1.e-6
                      * u.cm**-3)
    assert len(mean_density) == 9
    assert len(velocity) == len(mean_density) + 1
