from os.path import join as pathjoin
import pytest

import pandas as pd
import pandas.util.testing as pdt

from tardis.simulation import Simulation
from tardis.io.config_reader import Configuration

config_files = {
        'lte': 'tardis/plasma/tests/data/plasma_test_config_lte.yml',
        'nlte': 'tardis/plasma/tests/data/plasma_test_config_nlte.yml',
        }


class TestPlasmas():

    name = 'plasma_full/'

    @pytest.fixture(
            scope="class",
            params=config_files.items(),
            ids=config_files.keys()
            )
    def simulation(self, request, atomic_data_fname):
        name = request.param[0]
        config = Configuration.from_yaml(request.param[1])
        config['atom_data'] = atomic_data_fname
        simulation = Simulation.from_config(config)
        simulation.run()
        simulation._test_name = name
        return simulation

    def test_levels(self, simulation, tardis_ref_data):
        name = simulation._test_name
        new_levels = simulation.plasma.get_value('level_number_density')

        old_levels = tardis_ref_data[pathjoin(self.name, name, 'levels')]
        pdt.assert_almost_equal(
            new_levels, old_levels)

    def test_trads(self, simulation, tardis_ref_data):
        name = simulation._test_name

        new_t_rads = pd.Series(simulation.model.t_rad.to('K').value)

        old_t_rads = tardis_ref_data[pathjoin(self.name, name, 't_rad')]
        pdt.assert_almost_equal(
            new_t_rads, old_t_rads)
