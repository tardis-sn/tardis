import os
import warnings

import pandas as pd
import pytest
from pandas.util import testing as pdt
from numpy.testing import assert_almost_equal
from tardis.io.config_reader import Configuration
from tardis.simulation import Simulation

ionization = [
    {"ionization": "nebular"},
    {"ionization": "lte"},
]

excitation = [{"excitation": "lte"}, {"excitation": "dilute-lte"}]

radiative_rates_type = [
    {"radiative_rates_type": "detailed", "w_epsilon": 1.0e-10},
    {"radiative_rates_type": "detailed"},
    {"radiative_rates_type": "blackbody"},
    {"radiative_rates_type": "dilute-blackbody"},
]

line_interaction_type = [
    {"line_interaction_type": "scatter"},
    {"line_interaction_type": "macroatom"},
    {"line_interaction_type": "downbranch"},
]

disable_electron_scattering = [
    {"disable_electron_scattering": True},
    {"disable_electron_scattering": False},
]

nlte = [
    {"nlte": {"species": ["He I"], "coronal_approximation": True}},
    {"nlte": {"species": ["He I"], "classical_nebular": True}},
    {"nlte": {"species": ["He I"]}},
]

initial_t_inner = [{"initial_t_inner": "10000 K"}]

initial_t_rad = [{"initial_t_rad": "10000 K"}]

helium_treatment = [
    {"helium_treatment": "recomb-nlte"},
    {"helium_treatment": "recomb-nlte", "delta_treatment": 0.5},
]

config_list = (
    ionization
    + excitation
    + radiative_rates_type
    + line_interaction_type
    + disable_electron_scattering
    + nlte
    + initial_t_inner
    + initial_t_rad
    + helium_treatment
)


def idfn(fixture_value):
    """
    This function creates a string from a dictionary.
    We use it to obtain a readable name for the config fixture.
    """
    return str(
        "-".join(["{}:{}".format(k, v) for k, v in fixture_value.items()])
    )


class TestPlasma(object):

    general_properties = [
        "beta_rad",
        "g_electron",
        "selected_atoms",
        "number_density",
        "t_electrons",
        "w",
        "t_rad",
        "beta_electron",
    ]
    partiton_properties = ["level_boltzmann_factor", "partition_function"]
    atomic_properties = [
        "excitation_energy",
        "lines",
        "lines_lower_level_index",
        "lines_upper_level_index",
        "atomic_mass",
        "ionization_data",
        "nu",
        "wavelength_cm",
        "f_lu",
        "metastability",
    ]
    ion_population_properties = [
        "delta",
        "previous_electron_densities",
        "phi",
        "ion_number_density",
        "electron_densities",
    ]
    level_population_properties = ["level_number_density"]
    radiative_properties = [
        "stimulated_emission_factor",
        "previous_beta_sobolev",
        "tau_sobolevs",
        "beta_sobolev",
        "transition_probabilities",
    ]
    j_blues_properties = ["j_blues", "j_blues_norm_factor", "j_blue_estimator"]
    input_properties = ["volume", "r_inner"]
    helium_nlte_properties = ["helium_population", "helium_population_updated"]

    combined_properties = (
        general_properties
        + partiton_properties
        + atomic_properties
        + ion_population_properties
        + level_population_properties
        + radiative_properties
        + j_blues_properties
        + input_properties
        + helium_nlte_properties
    )

    scalars_properties = ["time_explosion", "link_t_rad_t_electron"]

    @pytest.fixture(scope="class")
    def chianti_he_db_fpath(self, tardis_ref_path):
        return os.path.abspath(
            os.path.join(tardis_ref_path, "atom_data", "chianti_He.h5")
        )

    @pytest.fixture(scope="class", params=config_list, ids=idfn)
    def config(self, request):
        config_path = os.path.join(
            "tardis", "plasma", "tests", "data", "plasma_base_test_config.yml"
        )
        config = Configuration.from_yaml(config_path)
        hash_string = ""
        for prop, value in request.param.items():
            hash_string = "_".join((hash_string, prop))
            if prop == "nlte":
                for nlte_prop, nlte_value in request.param[prop].items():
                    config.plasma.nlte[nlte_prop] = nlte_value
                    if nlte_prop != "species":
                        hash_string = "_".join((hash_string, nlte_prop))
            else:
                config.plasma[prop] = value
                hash_string = "_".join((hash_string, str(value)))
        hash_string = os.path.join("plasma_unittest", hash_string)
        setattr(config.plasma, "save_path", hash_string)
        return config

    @pytest.fixture(scope="class")
    def plasma(self, request, chianti_he_db_fpath, config, tardis_ref_data):
        config["atom_data"] = chianti_he_db_fpath
        sim = Simulation.from_config(config)
        if request.config.getoption("--generate-reference"):
            sim.plasma.to_hdf(tardis_ref_data, path=config.plasma.save_path)
            pytest.skip("Reference data saved at {0}".format(tardis_ref_data))
        return sim.plasma

    @pytest.mark.parametrize("attr", combined_properties)
    def test_plasma_properties(self, plasma, tardis_ref_data, config, attr):
        if hasattr(plasma, attr):
            actual = getattr(plasma, attr)
            if actual.ndim == 1:
                actual = pd.Series(actual)
            else:
                actual = pd.DataFrame(actual)
            key = os.path.join(config.plasma.save_path, "plasma", attr)
            expected = tardis_ref_data[key]
            pdt.assert_almost_equal(actual, expected)
        else:
            warnings.warn('Property "{}" not found'.format(attr))

    def test_levels(self, plasma, tardis_ref_data, config):
        actual = pd.DataFrame(plasma.levels)
        key = os.path.join(config.plasma.save_path, "plasma", "levels")
        expected = tardis_ref_data[key]
        pdt.assert_almost_equal(actual, expected)

    @pytest.mark.parametrize("attr", scalars_properties)
    def test_scalars_properties(self, plasma, tardis_ref_data, config, attr):
        actual = getattr(plasma, attr)
        if hasattr(actual, "cgs"):
            actual = actual.cgs.value
        key = os.path.join(config.plasma.save_path, "plasma", "scalars")
        expected = tardis_ref_data[key][attr]
        pdt.assert_almost_equal(actual, expected)

    def test_helium_treatment(self, plasma, tardis_ref_data, config):
        actual = plasma.helium_treatment
        key = os.path.join(config.plasma.save_path, "plasma", "scalars")
        expected = tardis_ref_data[key]["helium_treatment"]
        assert actual == expected

    def test_zeta_data(self, plasma, tardis_ref_data, config):
        if hasattr(plasma, "zeta_data"):
            actual = plasma.zeta_data
            key = os.path.join(config.plasma.save_path, "plasma", "zeta_data")
            expected = tardis_ref_data[key]
            assert_almost_equal(actual, expected.values)
