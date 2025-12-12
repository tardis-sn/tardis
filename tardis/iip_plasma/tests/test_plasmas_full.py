import os

import h5py
import numpy as np
import pytest
import yaml
from astropy import units as u

import tardis
from tardis.base import run_tardis

pytestmark = pytest.mark.skip("Skipping tests due to old format")


def data_path(fname):
    return os.path.join(
        tardis.__path__[0], "iip_plasma", "tests", "data", fname
    )


@pytest.fixture()
def plasma_compare_data_fname():
    return data_path("plasma_test_data.h5")


@pytest.fixture()
def plasma_compare_data(plasma_compare_data_fname):
    return h5py.File(plasma_compare_data_fname, "r")


class TestPlasmas:
    @classmethod
    @pytest.fixture(scope="class", autouse=True)
    def setup(self, request):
        atomic_dataset = request.config.getoption("--atomic-dataset")
        if not atomic_dataset:
            pytest.skip("--atomic-dataset was not specified")

        self.atom_data_filename = os.path.expanduser(
            os.path.expandvars(atomic_dataset)
        )
        assert os.path.exists(self.atom_data_filename), (
            f"{self.atom_data_filename} atomic datafiles does not seem to exist"
        )
        self.config_yaml = yaml.load(
            open("tardis/iip_plasma/tests/data/plasma_test_config_lte.yml"),
            Loader=yaml.CLoader,
        )
        self.config_yaml["atom_data"] = self.atom_data_filename
        self.lte_model = run_tardis(self.config_yaml)
        self.config_yaml = yaml.load(
            open("tardis/iip_plasma/tests/data/plasma_test_config_nlte.yml"),
            Loader=yaml.CLoader,
        )
        self.config_yaml["atom_data"] = self.atom_data_filename
        self.nlte_model = run_tardis(self.config_yaml)

    def test_lte_plasma(self, plasma_compare_data):
        old_plasma_t_rads = plasma_compare_data["test_lte1/t_rad"]
        old_plasma_levels = plasma_compare_data["test_lte1/levels"]

        new_plasma_t_rads = self.lte_model.t_rads / u.Unit("K")
        new_plasma_levels = (
            self.lte_model.plasma_array.get_value("level_number_density")
            .iloc[8]
            .iloc[1][10]
            .values
        )
        np.testing.assert_allclose(
            new_plasma_t_rads, old_plasma_t_rads, atol=100
        )
        np.testing.assert_allclose(
            new_plasma_levels, old_plasma_levels, rtol=0.1
        )

    def test_nlte_plasma(self, plasma_compare_data):
        old_plasma_t_rads = plasma_compare_data["test_nlte1/t_rad"]
        old_plasma_levels = plasma_compare_data["test_nlte1/levels"]
        new_plasma_t_rads = self.nlte_model.t_rads / u.Unit("K")
        new_plasma_levels = (
            self.nlte_model.plasma_array.get_value("level_number_density")
            .loc[2]
            .loc[1][10]
            .values
        )
        np.testing.assert_allclose(
            new_plasma_t_rads, old_plasma_t_rads, atol=150
        )
        np.testing.assert_allclose(
            new_plasma_levels, old_plasma_levels, rtol=0.1
        )
