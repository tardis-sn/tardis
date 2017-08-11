import os

import pandas as pd
import pytest
from pandas.util import testing as pdt
from numpy.testing import assert_almost_equal
from tardis.io.config_reader import Configuration
from tardis.simulation import Simulation


class BasePlasmaTest(object):
    #Class defining all common tests for different setups of Plasma
    #This can then be inherited for different Plasma setup
    @classmethod
    @pytest.fixture(scope="class", autouse=True)
    def setup(cls, atomic_data_fname, tardis_ref_path, reference_file_path, config_path):
        cls.config = Configuration.from_yaml(config_path)
        cls.config['atom_data'] = atomic_data_fname
        cls.sim = Simulation.from_config(cls.config)

        if pytest.config.getvalue("--generate-reference"):
            if os.path.exists(reference_file_path):
                pytest.skip(
                    'Reference data {0} does exist and tests will not '
                    'proceed generating new data'.format(reference_file_path))
            cls.sim.plasma.to_hdf(reference_file_path)
            pytest.skip("Reference data saved at {0}".format(
                reference_file_path))
        cls.plasma = cls.sim.plasma

    @pytest.fixture(scope="class")
    def reference_file_path(self):
        pass

    @pytest.fixture(scope="class")
    def config_path(self):
        pass

    #GENERAL PROPERTIES
    general_properties = ['beta_rad', 'g_electron', 'selected_atoms',
                          'number_density', 't_electrons', 'w', 't_rad']

    @pytest.mark.parametrize("attr", general_properties)
    def test_general_properties(self, reference_file_path, attr):
        actual = getattr(self.plasma, attr)
        expected = pd.read_hdf(reference_file_path,
                               os.path.join('plasma', attr))
        assert_almost_equal(actual, expected.values)

    #PARTITION PROPERTIES
    partiton_properties = ['level_boltzmann_factor', 'partition_function']

    @pytest.mark.parametrize("attr", partiton_properties)
    def test_partition_properties(self, reference_file_path, attr):
        actual = getattr(self.plasma, attr)
        expected = pd.read_hdf(reference_file_path,
                               os.path.join('plasma', attr))
        assert_almost_equal(actual, expected.values)

    #ATOMIC PROPERTIES
    atomic_properties = ['excitation_energy', 'lines', 'lines_lower_level_index',
                         'lines_upper_level_index', 'atomic_mass', 'ionization_data',
                         'nu', 'wavelength_cm', 'f_lu', 'metastability']

    @pytest.mark.parametrize("attr", atomic_properties)
    def test_atomic_properties(self, reference_file_path, attr):
        actual = getattr(self.plasma, attr)
        expected = pd.read_hdf(reference_file_path,
                               os.path.join('plasma', attr))
        assert_almost_equal(actual, expected.values)

    def test_levels(self, reference_file_path):
        actual = pd.DataFrame(self.plasma.levels)
        expected = pd.read_hdf(reference_file_path,
                               os.path.join('plasma', 'levels'))
        pdt.assert_almost_equal(actual, expected)

    #ION POPULATION PROPERTIES
    ion_population_properties = [
        'phi', 'ion_number_density', 'electron_densities']

    @pytest.mark.parametrize("attr", ion_population_properties)
    def test_ion_population_properties(self, reference_file_path, attr):
        actual = getattr(self.plasma, attr)
        expected = pd.read_hdf(reference_file_path,
                               os.path.join('plasma', attr))
        assert_almost_equal(actual, expected.values)

    #LEVEL POPULATION PROPERTIES
    level_population_properties = ['level_number_density']

    @pytest.mark.parametrize("attr", level_population_properties)
    def test_level_population_properties(self, reference_file_path, attr):
        actual = getattr(self.plasma, attr)
        expected = pd.read_hdf(reference_file_path,
                               os.path.join('plasma', attr))
        assert_almost_equal(actual, expected.values)

    #RADIATIVE PROPERTIES
    radiative_properties = ['stimulated_emission_factor',
                            'tau_sobolevs', 'beta_sobolev', 'transition_probabilities']

    @pytest.mark.parametrize("attr", radiative_properties)
    def test_radiative_properties(self, reference_file_path, attr):
        actual = getattr(self.plasma, attr)
        expected = pd.read_hdf(reference_file_path,
                               os.path.join('plasma', attr))
        assert_almost_equal(actual, expected.values)

    #J_BLUES PROPERTIES
    j_blues_properties = ['j_blues']

    @pytest.mark.parametrize("attr", j_blues_properties)
    def test_j_blues_properties(self, reference_file_path, attr):
        actual = getattr(self.plasma, attr)
        expected = pd.read_hdf(reference_file_path,
                               os.path.join('plasma', attr))
        assert_almost_equal(actual, expected.values)

    #SCALARS PROPERTIES
    scalars_properties = ['time_explosion', 'link_t_rad_t_electron']

    @pytest.mark.parametrize("attr", scalars_properties)
    def test_scalar_properties(self, reference_file_path, attr):
        actual = getattr(self.plasma, attr)
        if hasattr(actual, 'cgs'):
            actual = actual.cgs.value
        path = os.path.join('plasma', 'scalars')
        expected = pd.read_hdf(reference_file_path,
                               os.path.join('plasma', 'scalars'))[attr]
        assert_almost_equal(actual, expected)

    def test_helium_treatment(self, reference_file_path):
        actual = self.plasma.helium_treatment
        expected = pd.read_hdf(reference_file_path, os.path.join(
            'plasma', 'scalars'))['helium_treatment']
        assert actual == expected


class TestLTEPlasma(BasePlasmaTest):

    @pytest.fixture(scope="class")
    def reference_file_path(self, tardis_ref_path):
        return os.path.join(tardis_ref_path, 'plasma_reference', 'plasma_lte_reference.h5')

    @pytest.fixture(scope="class")
    def config_path(self):
        return os.path.join('tardis', 'plasma', 'tests', 'data', 'plasma_test_config_lte.yml')


class TestNLTEPlasma(BasePlasmaTest):

    @pytest.fixture(scope="class")
    def reference_file_path(self, tardis_ref_path):
        return os.path.join(tardis_ref_path, 'plasma_reference', 'plasma_nlte_reference.h5')

    @pytest.fixture(scope="class")
    def config_path(self):
        return os.path.join('tardis', 'plasma', 'tests', 'data', 'plasma_test_config_nlte.yml')

    #Additional Tests for NLTE Plasma
    #GENERAL PROPERTIES
    nlte_general_properties = ['beta_electron']

    @pytest.mark.parametrize("attr", nlte_general_properties)
    def test_nlte_general_properties(self, reference_file_path, attr):
        actual = getattr(self.plasma, attr)
        expected = pd.read_hdf(reference_file_path,
                               os.path.join('plasma', attr))
        assert_almost_equal(actual, expected.values)

    #ATOMIC PROPERTIES
    nlte_atomic_properties = ['zeta_data']

    @pytest.mark.parametrize("attr", nlte_atomic_properties)
    def test_nlte_atomic_properties(self, reference_file_path, attr):
        actual = getattr(self.plasma, attr)
        expected = pd.read_hdf(reference_file_path,
                               os.path.join('plasma', attr))
        assert_almost_equal(actual, expected.values)

    #ION POPULATION PROPERTIES
    nlte_ion_population_properties = ['delta', 'previous_electron_densities']

    @pytest.mark.parametrize("attr", nlte_ion_population_properties)
    def test_nlte_ion_population_properties(self, reference_file_path, attr):
        actual = getattr(self.plasma, attr)
        expected = pd.read_hdf(reference_file_path,
                               os.path.join('plasma', attr))
        assert_almost_equal(actual, expected.values)

    #RADIATIVE PROPERTIES
    nlte_radiative_properties = ['previous_beta_sobolev']

    @pytest.mark.parametrize("attr", nlte_radiative_properties)
    def test_nlte_radiative_properties(self, reference_file_path, attr):
        actual = getattr(self.plasma, attr)
        expected = pd.read_hdf(reference_file_path,
                               os.path.join('plasma', attr))
        assert_almost_equal(actual, expected.values)
