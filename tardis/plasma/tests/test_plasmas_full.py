import os
import pytest
import warnings
import copy

import pandas as pd
import pandas.util.testing as pdt

from tardis.simulation import Simulation
from tardis.io.config_reader import Configuration

config_files = {
    "lte": "tardis/plasma/tests/data/plasma_test_config_lte.yml",
    "nlte": "tardis/plasma/tests/data/plasma_test_config_nlte.yml",
}


class TestPlasmas:

    name = "plasma_full/"

    @pytest.fixture(scope="class")
    def refdata(self, tardis_ref_data):
        def get_ref_data(key):
            return tardis_ref_data[
                os.path.join(self.name, self._test_name, key)
            ]

        return get_ref_data

    @pytest.fixture(
        scope="class", params=config_files.items(), ids=config_files.keys()
    )
    def simulation(
        self, request, atomic_data_fname, generate_reference, tardis_ref_data
    ):
        name = request.param[0]
        config = Configuration.from_yaml(request.param[1])
        config["atom_data"] = atomic_data_fname
        simulation = Simulation.from_config(config)
        simulation.run()
        self._test_name = name

        if not generate_reference:
            return simulation
        else:
            simulation.plasma.hdf_properties = [
                "level_number_density",
            ]
            simulation.model.hdf_properties = ["t_radiative"]
            simulation.plasma.to_hdf(
                tardis_ref_data, self.name, self._test_name
            )
            simulation.model.to_hdf(tardis_ref_data, self.name, self._test_name)
            pytest.skip("Reference data was generated during this run.")
        return simulation

    def test_levels(self, simulation, refdata):
        new_levels = simulation.plasma.get_value("level_number_density")

        old_levels = refdata("level_number_density")
        pdt.assert_almost_equal(new_levels, old_levels)

    def test_trads(self, simulation, refdata):
        new_t_rads = pd.Series(simulation.model.t_rad.to("K").value)

        old_t_rads = refdata("t_radiative")
        pdt.assert_almost_equal(new_t_rads, old_t_rads)
