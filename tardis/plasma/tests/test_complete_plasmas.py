import os

import pandas as pd
import pytest
from pandas.util import testing as pdt
from numpy.testing import assert_almost_equal
from tardis.io.config_reader import Configuration
from tardis.simulation import Simulation

ionization = [
    {'ionization': 'nebular'},
    {'ionization': 'lte'},
]

excitation = [
    {'excitation': 'lte'},
    {'excitation': 'dilute-lte'}
]

radiative_rates_type = [
    {'radiative_rates_type': 'detailed', 'w_epsilon': 1.0e-10},
    {'radiative_rates_type': 'detailed'},
    {'radiative_rates_type': 'blackbody'},
    {'radiative_rates_type': 'dilute-blackbody'}
]

line_interaction_type = [
    {'line_interaction_type': 'scatter'},
    {'line_interaction_type': 'macroatom'},
    {'line_interaction_type': 'downbranch'}
]

disable_electron_scattering = [
    {'disable_electron_scattering': True},
    {'disable_electron_scattering': False}
]

nlte = [
    {'nlte':  {'species': ['He I'], 'coronal_approximation': True}},
    {'nlte':  {'species': ['He I'], 'classical_nebular': True}},
    {'nlte':  {'species': ['He I']}}
]

initial_t_inner = [
    {'initial_t_inner': '10000 K'}
]

initial_t_rad = [
    {'initial_t_rad': '10000 K'}
]

helium_treatment = [
    {'helium_treatment': 'recomb-nlte'},
    {'helium_treatment': 'recomb-nlte', 'delta_treatment': 0.5}
]

config_list = (
        ionization + excitation + radiative_rates_type +
        line_interaction_type + disable_electron_scattering + nlte +
        initial_t_inner + initial_t_rad + helium_treatment)


def idfn(fixture_value):
    return str('-'.join([
        '{}:{}'.format(k, v) for k, v in fixture_value.items()]))


class TestPlasma(object):

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
    helium_nlte_properties = ['helium_population', 'helium_population_updated']

    combined_properties = (
            general_properties + partiton_properties +
            atomic_properties + ion_population_properties +
            level_population_properties + radiative_properties +
            j_blues_properties + input_properties + helium_nlte_properties)

    scalars_properties = ['time_explosion', 'link_t_rad_t_electron']

    @pytest.fixture(scope="class")
    def chianti_he_db_fpath(self, tardis_ref_path):
        return os.path.abspath(os.path.join(
            tardis_ref_path, 'atom_data', 'chianti_He.h5'))

    @pytest.fixture(scope="class")
    def reference_fpath(self, tardis_ref_path):
        path = os.path.join(
            tardis_ref_path, 'plasma_reference', 'reference_data.h5')
        if pytest.config.getvalue("--generate-reference"):
            if os.path.exists(path):
                pytest.skip('Reference data {0} does exist and tests will not '
                            'proceed generating new data'.format(path))
        return path

    @pytest.fixture(
            scope="class",
            params=config_list,
            ids=idfn
            )
    def config(self, request):
        config_path = os.path.join(
            'tardis', 'plasma', 'tests', 'data', 'plasma_base_test_config.yml')
        config = Configuration.from_yaml(config_path)
        hash_string = ''
        for prop, value in request.param.items():
            hash_string = '_'.join((hash_string, prop))
            if prop == 'nlte':
                for nlte_prop, nlte_value in request.param[prop].items():
                    config.plasma.nlte[nlte_prop] = nlte_value
                    if nlte_prop != 'species':
                        hash_string = '_'.join((hash_string, nlte_prop))
            else:
                config.plasma[prop] = value
                hash_string = '_'.join((hash_string, str(value)))
        setattr(config.plasma, 'save_path', hash_string)
        return config

    @pytest.yield_fixture(scope="class")
    def reference(self, reference_fpath):
        with pd.HDFStore(reference_fpath) as hdf_file:
            yield hdf_file

    @pytest.fixture(scope="class")
    def plasma(self, chianti_he_db_fpath, config, reference_fpath, reference):
        config['atom_data'] = chianti_he_db_fpath
        sim = Simulation.from_config(config)
        if pytest.config.getvalue("--generate-reference"):
            sim.plasma.to_hdf(reference_fpath, path=config.plasma.save_path)
            pytest.skip("Reference data saved at {0}".format(reference_fpath))
        return sim.plasma


    @pytest.mark.parametrize("attr", combined_properties)
    def test_plasma_properties(self, plasma, reference, config, attr):
        if hasattr(plasma, attr):
            actual = getattr(plasma, attr)
            if actual.ndim == 1:
                actual = pd.Series(actual)
            else:
                actual = pd.DataFrame(actual)
            expected = reference.select(os.path.join(
                config.plasma.save_path, 'plasma', attr))
            pdt.assert_almost_equal(actual, expected)

    def test_levels(self, plasma, reference, config):
        actual = pd.DataFrame(plasma.levels)
        expected = reference.select(os.path.join(
            config.plasma.save_path, 'plasma', 'levels'))
        pdt.assert_almost_equal(actual, expected)

    @pytest.mark.parametrize("attr", scalars_properties)
    def test_scalars_properties(self, plasma, reference, config, attr):
        actual = getattr(plasma, attr)
        if hasattr(actual, 'cgs'):
            actual = actual.cgs.value
        expected = reference.select(os.path.join(
            config.plasma.save_path, 'plasma', 'scalars'))[attr]
        pdt.assert_almost_equal(actual, expected)

    def test_helium_treatment(self, plasma, reference, config):
        actual = plasma.helium_treatment
        expected = reference.select(os.path.join(
            config.plasma.save_path, 'plasma', 'scalars'))['helium_treatment']
        assert actual == expected

    def test_zeta_data(self, plasma, reference, config):
        if hasattr(plasma, 'zeta_data'):
            actual = plasma.zeta_data
            expected = reference.select(os.path.join(
                config.plasma.save_path, 'plasma', 'zeta_data'))
            assert_almost_equal(actual, expected.values)
