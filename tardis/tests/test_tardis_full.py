import pytest
import os
import yaml
from tardis import io, model, simulation
from numpy.testing import assert_almost_equal, assert_array_almost_equal

@pytest.mark.pure_kurucz
class TestSimpleRun1():
    """
    Very simple run with excitation and ionization set to
    """
    @classmethod
    @pytest.fixture(scope = "class", autouse = True)
    def setup(self, pure_kurucz_filename):
        self.atom_data_filename = os.path.expanduser(pure_kurucz_filename)
        assert os.path.exists(self.atom_data_filename)
        self.config_yaml = yaml.load(open('tardis/io/tests/data/tardis_configv1_verysimple.yml'))
        self.config_yaml['atom_data'] = self.atom_data_filename

        self.config = io.config_reader.TARDISConfiguration.from_config_dict(self.config_yaml)
        assert self.config.atom_data.uuid1 == 'ca2684c53cc511e39f0ec8bcc8a04795', 'requiring specific atom dataset'
        assert self.config.atom_data.md5 == '4cccdc9b4c09faf3e348a2c710fb1715', 'requiring specific atom dataset'
        self.model = model.Radial1DModel(self.config)
        simulation.run_radial1d(self.model)


    def test_structure(self):
        assert_array_almost_equal(self.config.v_inner, 0)