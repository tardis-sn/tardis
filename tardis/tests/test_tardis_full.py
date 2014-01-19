import pytest
import os
import yaml
from tardis import io


@pytest.mark.pure_kurucz
class TestPaper1Run():

    @classmethod
    @pytest.fixture(scope = "class", autouse = True)
    def setup(self, pure_kurucz_filename):
        self.atom_data_filename = os.path.expanduser(pure_kurucz_filename)
        assert os.path.exists(self.atom_data_filename)
        self.config_yaml = yaml.load(open('tardis/io/tests/data/paper1_tardis_configv1.yml'))
        self.config_yaml['atom_data'] = self.atom_data_filename

        self.config = io.config_reader.TARDISConfiguration.from_config

    def test_first(self):
        print self.atom_data_filename