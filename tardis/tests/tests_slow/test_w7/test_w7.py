import os
import yaml
import numpy as np
import pytest
from numpy.testing import assert_allclose
from astropy import units as u

from tardis.simulation.base import Simulation
from tardis.model import Radial1DModel
from tardis.io.config_reader import Configuration


def data_path(fname):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fname)


@pytest.mark.skipif(not pytest.config.getvalue("atomic-dataset"),
                    reason='--atomic_database was not specified')
@pytest.mark.skipif(not pytest.config.getvalue("with-slow"),
                    reason="slow tests can only be run using --with-slow")
class SlowTestW7(object):
    """
    Slow integration test for Stratified W7 setup.

    Assumed two compressed binaries (.npz) are placed in this directory:

    * baseline_ndarrays.npz              | * baseline_quantities.npz
    Contents (all (.npy)):               | Contents (all (.npy)):
        * last_interaction_type          |     * t_rads
        * last_line_interaction_out_id   |     * luminosity_inner
        * last_line_interaction_in_id    |     * montecarlo_luminosity
        * j_estimators                   |     * montecarlo_virtual_luminousity
        * j_blue_estimators              |     * time_of_simulation
        * last_line_interaction_shell_id |     * montecarlo_nu
        * nubar_estimators               |     * last_line_interaction_angstrom
        * ws                             |     * j_blues_norm_factor

    Also assumed `kurucz_cd23_chianti_H_He.h5` file exists in `tmp` directory.
    """

    @classmethod
    @pytest.fixture(scope="class", autouse=True)
    def setup(self):
        """
        This method does initial setup of creating configuration and performing
        a single run of integration test.
        """
        self.config_file = data_path("config_w7.yml")
        self.abundances = data_path("abundancies_w7.dat")
        self.densities = data_path("densities_w7.dat")

        # First we check whether all files exist at desired paths.
        assert os.path.exists(self.config_file), \
            '%s does not exist' % self.config_file

        assert os.path.exists(self.abundances), \
            '%s does not exist' % self.abundances

        assert os.path.exists(self.densities), \
            '%s does not exist' % self.densities

        assert os.path.exists('/tmp/kurucz_cd23_chianti_H_He.h5'), \
            'kurucz_cd23_chianti_H_He.h5 atom data file does not exist'

        # The available config file doesn't have file paths of atom data file,
        # densities and abundances profile files as desired. We form dictionary
        # from the config file and override those parameters by putting file
        # paths of these three files at proper places.
        config_yaml = yaml.load(open(self.config_file))
        config_yaml['atom_data'] = '/tmp/kurucz_cd23_chianti_H_He.h5'
        config_yaml['model']['abundances']['filename'] = self.abundances
        config_yaml['model']['structure']['filename'] = self.densities

        # The config hence obtained will be having appropriate file paths.
        tardis_config = Configuration.from_config_dict(config_yaml)

        # We now do a run with prepared config and get radial1d model.
        self.obtained_radial1d_model = Radial1DModel(tardis_config)
        simulation = Simulation(tardis_config)
        simulation.legacy_run_simulation(self.obtained_w7_radial1d_model)

        # The baseline data against which assertions are to be made is ingested
        # from already available compressed binaries (.npz). These will return
        # dictionaries of numpy.ndarrays for performing assertions.
        self.expected_ndarrays = np.load(data_path("baseline_ndarrays.npz"))
        self.expected_quantities = np.load(data_path("baseline_quantities.npz"))
