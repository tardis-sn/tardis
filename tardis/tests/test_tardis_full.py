import pytest
import os
import yaml
from tardis import io, model, simulation
from tardis.io.config_reader import TARDISConfiguration
from numpy.testing import assert_array_almost_equal

@pytest.mark.skipif(not pytest.config.getvalue("atomic-dataset"),
                    reason='--atomic_database was not specified')
class TestSimpleRun():
    """
    Very simple run
    """

    @classmethod
    @pytest.fixture(scope="class", autouse=True)
    def setup(self):
        self.atom_data_filename = pytest.config.getvalue('atomic-dataset')
        assert os.path.exists(self.atom_data_filename)
        self.config_yaml = yaml.load(open('tardis/io/tests/data/tardis_configv1_verysimple.yml'))
        self.config_yaml['atom_data'] = self.atom_data_filename

        self.config = TARDISConfiguration.from_config_dict(self.config_yaml)
        self.model = model.Radial1DModel(self.config)
        simulation.run_radial1d(self.model)


    def test_structure(self):
        pass