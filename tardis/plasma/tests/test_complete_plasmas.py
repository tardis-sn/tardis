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

    @pytest.fixture(scope="class")
    def plasma(self, chianti_he_db_fpath, tardis_ref_path, data):
        config = data['config']
        config['atom_data'] = chianti_he_db_fpath
        sim = Simulation.from_config(config)

        if pytest.config.getvalue("--generate-reference"):
            if os.path.exists(data['reference_file_path']):
                pytest.skip(
                    'Reference data {0} does exist and tests will not '
                    'proceed generating new data'.format(data['reference_file_path']))
            sim.plasma.to_hdf(data['reference_file_path'])
            pytest.skip("Reference data saved at {0}".format(
                data['reference_file_path']))
        return sim.plasma

    @pytest.fixture(scope="class")
    def chianti_he_db_fpath(self):
        return os.path.abspath(os.path.join('tardis', 'tests', 'data', 'chianti_he_db.h5'))

    @pytest.fixture(scope="class")
    def data(self):
        pass

    @pytest.yield_fixture()
    def reference(self, data):
        with pd.HDFStore(data['reference_file_path']) as hdf_file:
            yield hdf_file

    general_properties = ['beta_rad', 'g_electron', 'selected_atoms',
                          'number_density', 't_electrons', 'w', 't_rad']
    partiton_properties = ['level_boltzmann_factor', 'partition_function']
    atomic_properties = ['excitation_energy', 'lines', 'lines_lower_level_index',
                         'lines_upper_level_index', 'atomic_mass', 'ionization_data',
                         'nu', 'wavelength_cm', 'f_lu', 'metastability']
    ion_population_properties = [
        'phi', 'ion_number_density', 'electron_densities']
    level_population_properties = ['level_number_density']
    radiative_properties = ['stimulated_emission_factor',
                            'tau_sobolevs', 'beta_sobolev', 'transition_probabilities']
    j_blues_properties = ['j_blues']

    combined_properties = general_properties + partiton_properties + atomic_properties + ion_population_properties + \
        level_population_properties + radiative_properties + \
        j_blues_properties

    @pytest.mark.parametrize("attr", combined_properties)
    def test_plasma_properties(self, plasma, reference, attr):
        actual = getattr(plasma, attr)
        if actual.ndim == 1:
            actual = pd.Series(actual)
        else:
            actual = pd.DataFrame(actual)
        expected = reference.select(os.path.join('plasma', attr))
        pdt.assert_almost_equal(actual, expected)

    def test_levels(self, plasma, reference):
        actual = pd.DataFrame(plasma.levels)
        expected = reference.select(os.path.join('plasma', 'levels'))
        pdt.assert_almost_equal(actual, expected)

    scalars_properties = ['time_explosion', 'link_t_rad_t_electron']

    @pytest.mark.parametrize("attr", scalars_properties)
    def test_scalars_properties(self, plasma, reference, attr):
        actual = getattr(plasma, attr)
        if hasattr(actual, 'cgs'):
            actual = actual.cgs.value
        expected = reference.select(os.path.join(
            'plasma', 'scalars'))[attr]
        pdt.assert_almost_equal(actual, expected)

    def test_helium_treatment(self, plasma, reference):
        actual = plasma.helium_treatment
        expected = reference.select(os.path.join(
            'plasma', 'scalars'))['helium_treatment']
        assert actual == expected


class TestLTEPlasma(BasePlasmaTest):

    @pytest.fixture(scope="class")
    def data(self, tardis_ref_path):
        plasma_data = {}
        plasma_data['reference_file_path'] = os.path.join(
            tardis_ref_path, 'plasma_reference', 'plasma_lte_reference.h5')
        config_path = os.path.join(
            'tardis', 'plasma', 'tests', 'data', 'plasma_test_config_lte.yml')
        plasma_data['config'] = Configuration.from_yaml(config_path)
        return plasma_data


class TestNLTEPlasma(BasePlasmaTest):

    @pytest.fixture(scope="class")
    def data(self, tardis_ref_path):
        plasma_data = {}
        plasma_data['reference_file_path'] = os.path.join(
            tardis_ref_path, 'plasma_reference', 'plasma_nlte_reference.h5')
        config_path = os.path.join(
            'tardis', 'plasma', 'tests', 'data', 'plasma_test_config_nlte.yml')
        plasma_data['config'] = Configuration.from_yaml(config_path)
        return plasma_data

    #Additional Tests for NLTE Plasma

    nlte_general_properties = ['beta_electron']
    nlte_ion_population_properties = ['delta', 'previous_electron_densities']
    nlte_radiative_properties = ['previous_beta_sobolev']
    nlte_properties = nlte_general_properties + \
        nlte_ion_population_properties + nlte_radiative_properties

    @pytest.mark.parametrize("attr", nlte_properties)
    def test_nlte_properties(self, plasma, reference, attr):
        actual = getattr(plasma, attr)
        if actual.ndim == 1:
            actual = pd.Series(actual)
        else:
            actual = pd.DataFrame(actual)
        expected = reference.select(os.path.join('plasma', attr))
        pdt.assert_almost_equal(actual, expected)

    def test_zeta_data(self, plasma, reference):
        actual = getattr(plasma, 'zeta_data')
        expected = reference.select(os.path.join('plasma', 'zeta_data'))
        assert_almost_equal(actual, expected.values)


class TestPlasmaSetupIII(TestNLTEPlasma):

    @pytest.fixture(scope="class")
    def data(self, tardis_ref_path):
        plasma_data = {}
        plasma_data['reference_file_path'] = os.path.join(
            tardis_ref_path, 'plasma_reference', 'plasma_setup_III_reference.h5')
        plasma_data['config_path'] = os.path.join(
            'tardis', 'plasma', 'tests', 'data', 'plasma_test_config_nlte.yml')
        config = Configuration.from_yaml(plasma_data['config_path'])
        config.plasma.radiative_rates_type = 'detailed'
        config.plasma.nlte.classical_nebular = True
        #config.plasma.helium_treatment = 'recomb-nlte'
        plasma_data['config'] = config
        return plasma_data

    j_blues_detailed_properties = ['j_blues_norm_factor', 'j_blue_estimator']
    additional_properties = j_blues_detailed_properties + ['volume', 'r_inner']

    @pytest.mark.parametrize("attr", additional_properties)
    def test_j_blues_detailed_properties(self, plasma, reference,  attr):
        actual = getattr(plasma, attr)
        if actual.ndim == 1:
            actual = pd.Series(actual)
        else:
            actual = pd.DataFrame(actual)
        expected = reference.select(os.path.join('plasma', attr))
        pdt.assert_almost_equal(actual, expected)


class TestPlasmaSetupIV(TestNLTEPlasma):

    @pytest.fixture(scope="class")
    def data(self, tardis_ref_path):
        plasma_data = {}
        plasma_data['reference_file_path'] = os.path.join(
            tardis_ref_path, 'plasma_reference', 'plasma_setup_IV_reference.h5')
        plasma_data['config_path'] = os.path.join(
            'tardis', 'plasma', 'tests', 'data', 'plasma_test_config_nlte.yml')
        config = Configuration.from_yaml(plasma_data['config_path'])
        config.plasma.radiative_rates_type = 'blackbody'
        config.plasma.nlte.coronal_approximation = True
        #config.plasma.helium_treatment = 'numerical-nlte'
        plasma_data['config'] = config
        return plasma_data
