import warnings
from pathlib import Path

import pandas as pd
import pytest
from numpy import testing as npt
from pandas import testing as pdt

from tardis.io.configuration.config_reader import Configuration
from tardis.simulation import Simulation
from tardis.tests.fixtures.regression_data import RegressionData


import unittest
import numpy as np
from tardis.plasma.properties.ion_population import get_zeta_values  # Adjust the import if necessary

class TestZetaValues(unittest.TestCase):
    def setUp(self):
        """
        This will run before every test case. Use it to prepare test data.
        """
        # Prepare mock zeta_data (this is a mock, you can adjust it according to your data)
        # Example: Columns are t_rads, and rows are ion indices
        data = {
            0.1: [0.5, 0.6, 0.7],
            0.2: [0.4, 0.5, 0.6],
            0.3: [0.3, 0.4, 0.5],
            0.4: [0.2, 0.3, 0.4]
        }
        self.zeta_data = pd.DataFrame(data, index=[1, 2, 3]) # Ion indices as row labels
        self.ion_index = 2  

    def test_get_zeta_values_with_active_shells(self):
        """
        Test the zeta_values function with active shells only.
        """
        active_shells = [True, False, True, True] 
        t_rad = np.array([0.1, 0.3, 0.4])  

        zeta_values = get_zeta_values(self.zeta_data, self.ion_index, t_rad, active_shells=active_shells)

        self.assertEqual(len(zeta_values), 3)  
        self.assertTrue(np.all(np.isfinite(zeta_values))) 

    def test_get_zeta_values_without_active_shells(self):
        """
        Test the zeta_values function without filtering active shells (i.e., using all shells).
        """
        t_rad = np.array([0.1, 0.2, 0.3, 0.4]) 

        zeta_values = get_zeta_values(self.zeta_data, self.ion_index, t_rad)

        self.assertEqual(len(zeta_values), 4) 
        self.assertTrue(np.all(np.isfinite(zeta_values)))  

    def test_get_zeta_values_with_out_of_range_t_rad(self):
        """
        Test the zeta_values function with t_rads outside the interpolation range.
        """
        t_rad = np.array([0.5]) 

        zeta_values = get_zeta_values(self.zeta_data, self.ion_index, t_rad)

      
        self.assertEqual(zeta_values[0], 1.0) 


PLASMA_CONFIG_FPATH = (
    Path("tardis") / "plasma" / "tests" / "data" / "plasma_base_test_config.yml"
)


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

CONFIG_LIST = (
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
    return "-".join([f"{k}:{v}" for k, v in fixture_value.items()])


class TestPlasma:
    regression_data = None

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
    def chianti_he_db_fpath(self, tardis_regression_path):
        return (
            tardis_regression_path / "atom_data" / "chianti_He.h5"
        ).absolute()

    @pytest.fixture(scope="class", params=CONFIG_LIST, ids=idfn)
    def config(self, request):
        config = Configuration.from_yaml(PLASMA_CONFIG_FPATH)
        hash_string = ""
        for prop, value in request.param.items():
            hash_string = f"{hash_string}_{prop}"
            if prop == "nlte":
                for nlte_prop, nlte_value in request.param[prop].items():
                    config.plasma.nlte[nlte_prop] = nlte_value
                    if nlte_prop != "species":
                        hash_string = f"{hash_string}_{nlte_prop}"
            else:
                config.plasma[prop] = value
                hash_string = "_".join((hash_string, str(value)))
        hash_string = f"plasma_unittest{hash_string}"
        config.plasma.save_path = hash_string
        request.cls.regression_data = RegressionData(request)
        request.cls.regression_data.fname = f"{hash_string}.h5"
        return config

    @pytest.fixture(scope="class")
    def plasma(
        self,
        chianti_he_db_fpath,
        config,
    ):
        config["atom_data"] = str(chianti_he_db_fpath)
        sim = Simulation.from_config(config)
        data = self.regression_data.sync_hdf_store(
            sim.plasma, update_fname=False
        )
        yield sim.plasma
        data.close()

    @pytest.mark.parametrize("attr", combined_properties)
    def test_plasma_properties(self, plasma, attr):
        key = f"plasma/{attr}"
        try:
            expected = pd.read_hdf(self.regression_data.fpath, key)
        except KeyError:
            pytest.skip(f"Key {key} not found in regression data")

        if hasattr(plasma, attr):
            actual = getattr(plasma, attr)
            if attr == "selected_atoms":
                npt.assert_allclose(actual.values, expected.values)
            elif actual.ndim == 1:
                actual = pd.Series(actual)
                pdt.assert_series_equal(actual, expected)
            else:
                actual = pd.DataFrame(actual)
                pdt.assert_frame_equal(actual, expected)
        else:
            warnings.warn(f'Property "{attr}" not found')

    def test_levels(self, plasma):
        actual = pd.DataFrame(plasma.levels)
        key = "plasma/levels"
        expected = pd.read_hdf(self.regression_data.fpath, key)
        pdt.assert_frame_equal(actual, expected)

    @pytest.mark.parametrize("attr", scalars_properties)
    def test_scalars_properties(self, plasma, attr):
        actual = getattr(plasma, attr)
        if hasattr(actual, "cgs"):
            actual = actual.cgs.value
        key = "plasma/scalars"
        expected = pd.read_hdf(self.regression_data.fpath, key)[attr]
        npt.assert_equal(actual, expected)

    def test_helium_treatment(self, plasma):
        actual = plasma.helium_treatment
        key = "plasma/scalars"
        expected = pd.read_hdf(self.regression_data.fpath, key)[
            "helium_treatment"
        ]
        assert actual == expected

    def test_zeta_data(self, plasma):
        if hasattr(plasma, "zeta_data"):
            actual = plasma.zeta_data
            key = "plasma/zeta_data"
            expected = pd.read_hdf(self.regression_data.fpath, key)
            pdt.assert_frame_equal(actual, expected)

