import os

import pandas as pd
import pytest
from pandas.util import testing as pdt
from numpy.testing import assert_almost_equal
from tardis.io.config_reader import Configuration
from tardis.simulation import Simulation


setupI = {
    'ionization': 'lte',
    'excitation': 'lte',
    'radiative_rates_type': 'dilute-blackbody',
    'line_interaction_type': 'scatter',
    'reference_file_path': 'plasma_setupI_reference.h5'
}

setupII = {
    'ionization': 'nebular',
    'excitation': 'dilute-lte',
    'radiative_rates_type': 'blackbody',
    'line_interaction_type': 'macroatom',
    'nlte': {'species': ['He I'], 'classical_nebular': True},
    'initial_t_inner': '10000 K',
    'initial_t_rad': '10000 K',
    'disable_electron_scattering': False,
    'w_epsilon': '1.0e-10',
    'reference_file_path': 'plasma_setupII_reference.h5'
}

setupIII = {
    'ionization': 'nebular',
    'excitation': 'dilute-lte',
    'radiative_rates_type': 'detailed',
    'line_interaction_type': 'downbranch',
    'disable_electron_scattering': True,
    'nlte': {'species': ['He I'], 'coronal_approximation': True},
    'reference_file_path': 'plasma_setupIII_reference.h5'
}


class TestPlasma(object):

    @pytest.fixture(scope="class")
    def plasma(self, chianti_he_db_fpath, config):
        config['atom_data'] = chianti_he_db_fpath
        sim = Simulation.from_config(config)
        if pytest.config.getvalue("--generate-reference"):
            if os.path.exists(config.reference_file_path):
                pytest.skip(
                    'Reference data {0} does exist and tests will not '
                    'proceed generating new data'.format(config.reference_file_path))
            sim.plasma.to_hdf(config.reference_file_path)
            pytest.skip("Reference data saved at {0}".format(
                config.reference_file_path))
        return sim.plasma

    @pytest.fixture(scope="class")
    def chianti_he_db_fpath(self):
        return os.path.abspath(os.path.join('tardis', 'tests', 'data', 'chianti_he_db.h5'))

    @pytest.fixture(scope="class", params=[setupI, setupII, setupIII])
    def config(self, request, tardis_ref_path):
        config_path = os.path.join(
            'tardis', 'plasma', 'tests', 'data', 'plasma_base_test_config.yml')
        config = Configuration.from_yaml(config_path)
        for prop, value in request.param.items():
            if prop == 'reference_file_path':
                setattr(config, prop, os.path.join(
                    tardis_ref_path, 'plasma_reference', value))
            elif prop == 'nlte':
                for nlte_prop, nlte_value in request.param[prop].items():
                    setattr(config.plasma.nlte, nlte_prop, nlte_value)
            else:
                setattr(config.plasma, prop, value)
        return config

    @pytest.yield_fixture()
    def reference(self, config):
        with pd.HDFStore(config.reference_file_path) as hdf_file:
            yield hdf_file

    general_properties = ['beta_rad', 'g_electron', 'selected_atoms',
                          'number_density', 't_electrons', 'w', 't_rad', 'beta_electron']
    partiton_properties = ['level_boltzmann_factor', 'partition_function']
    atomic_properties = ['excitation_energy', 'lines', 'lines_lower_level_index',
                         'lines_upper_level_index', 'atomic_mass', 'ionization_data',
                         'nu', 'wavelength_cm', 'f_lu', 'metastability']
    ion_population_properties = ['delta', 'previous_electron_densities',
                                 'phi', 'ion_number_density', 'electron_densities']
    level_population_properties = ['level_number_density']
    radiative_properties = ['stimulated_emission_factor', 'previous_beta_sobolev',
                            'tau_sobolevs', 'beta_sobolev', 'transition_probabilities']
    j_blues_properties = ['j_blues', 'j_blues_norm_factor', 'j_blue_estimator']
    input_properties = ['volume', 'r_inner']

    combined_properties = general_properties + partiton_properties + atomic_properties + \
        ion_population_properties + level_population_properties + radiative_properties + \
        j_blues_properties + input_properties

    @pytest.mark.parametrize("attr", combined_properties)
    def test_plasma_properties(self, plasma, reference, attr):
        if hasattr(plasma, attr):
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

    def test_zeta_data(self, plasma, reference):
        if hasattr(plasma, 'zeta_data'):
            actual = plasma.zeta_data
            expected = reference.select(os.path.join('plasma', 'zeta_data'))
            assert_almost_equal(actual, expected.values)
