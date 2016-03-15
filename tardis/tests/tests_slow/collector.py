import os
import yaml
import numpy as np
from tardis import run_tardis


def data_path(model_dir, filename):
    """Takes in strings and returns a decorated file path.

    Args:
        model_dir: (~str)
            Particular directory name in `tests_slow` directory containing files
            for a specific setup.
        fname: (~str)
            Name of file inside model_dir.

    Returns:
        (~str) $TARDIS_PATH/tests/tests_slow/model_dir/fname
    """
    return os.path.join(os.path.join(os.path.dirname(os.path.realpath(__file__)), model_dir, filename))


class Collector:
    """Collect the benchmark data by performing a single run of provided setup. This
    benchmark data will be kept for future assertions. This class is only to be used
    by trusted version of TARDIS which is certain to provide correct output Radial1D
    model from the run.

    Args:
        model_dir: (~str)
            Particular directory name in `tests_slow` directory containing files
            for a specific setup.
        config: (~str or ~dict)
            Name of yml config file inside model_dir, or a dict containing config.
        abundances: (~str)
            Name of the data file containing abundances profile.
        densities: (~str)
            Name of the data file containing densities profile.
        radial1d_model: (~tardis.model.Radial1DModel)
            Initialized to None first - will contain the information after performing
            the run.
    """

    def __init__(self, model_dir, config, abundances=None, densities=None, radial1d_model=None):
        self.model_dir = model_dir
        self.config = config
        self.abundances = abundances
        self.densities = densities
        self.radial1d_model = radial1d_model

    @staticmethod
    def save_kurucz_atom_data_file():
        """Download `kurucz_cd23_chianti_H_He.h5` and save it in /tmp directory."""

        # To keep the tardis working repo clean from a large atom data file, we store
        # this file in the 'tmp' directory, and notify the user about its existence
        if not os.path.exists(os.path.abspath('/tmp/kurucz_cd23_chianti_H_He.h5')):
            os.system('wget http://www.mpa-garching.mpg.de/~michi/tardis/data/kurucz_cd23_chianti_H_He.zip')
            os.system('unzip kurucz_cd23_chianti_H_He.zip')
            os.system('mv kurucz_cd23_chianti_H_He.h5 /tmp/')
            os.system('rm kurucz_cd23_chianti_H_He.zip')
            print('Successfully obtained atom data file and store in /tmp directory.')
        else:
            print('The atom data file already exists in /tmp directory.')

    def perform_run(self):
        """Performs a single run of TARDIS with the specified model, stores the results in
        self.radial1d_model.
        """
        # Path to the yaml config files and atom data files
        yml_config_filepath = data_path(self.model_dir, self.config)
        self.save_kurucz_atom_data_file()

        # Obtain a dictionary out of yaml config filepath and override values of the
        # path to atom data, densities and/or abundances data files if provided
        config_dict = yaml.load(open(yml_config_filepath))
        config_dict['atom_data'] = os.path.abspath('/tmp/kurucz_cd23_chianti_H_He.h5')
        if self.abundances is not None:
            config_dict['model']['abundances']['filename'] = data_path(self.model_dir, self.abundances)

        if self.densities is not None:
            config_dict['model']['structure']['filename'] = data_path(self.model_dir, self.densities)

        # Prepare a standard model from a 'believed' stable build of TARDIS
        self.radial1d_model = run_tardis(config_dict)

    def collect_ndarrays_from_radial1d_model(self, filename='expected_ndarrays.npz'):
        """Save all the members of radial1d_model which are numpy.ndarrays into a single
        file in compressed binary (`.npz`) format.

        Args:
            filename : (~str)
                File name of `.npz` file.
        """
        np.savez_compressed(
            data_path(self.model_dir, filename),
            last_interaction_type=self.radial1d_model.last_interaction_type,
            last_line_interaction_out_id=self.radial1d_model.last_line_interaction_out_id,
            last_line_interaction_in_id=self.radial1d_model.last_line_interaction_in_id,
            j_estimators=self.radial1d_model.j_estimators,
            j_blue_estimators=self.radial1d_model.j_blue_estimators,
            last_line_interaction_shell_id=self.radial1d_model.last_line_interaction_shell_id,
            nubar_estimators=self.radial1d_model.nubar_estimators,
            ws=self.radial1d_model.ws)

    def collect_astropy_quantities_from_radial1d_model(self, filename='expected_quantities.npz'):
        """Save all the members of radial1d_model which are numpy.ndarrays into a single
        file in compressed binary (`.npz`) format. Only values are stored and while just
        before performing assertions, the model under testing is converted into the units
        same as what benchmark data had earlier.

        Args:
            filename : (~str)
                File name of `.npz` file.
        """
        np.savez_compressed(
            data_path(self.model_dir, filename),
            t_rads=self.radial1d_model.t_rads.value,
            luminosity_inner=self.radial1d_model.luminosity_inner.value,
            montecarlo_luminosity=self.radial1d_model.montecarlo_luminosity.value,
            montecarlo_virtual_luminosity=self.radial1d_model.montecarlo_virtual_luminosity.value,
            time_of_simulation=self.radial1d_model.time_of_simulation.value,
            montecarlo_nu=self.radial1d_model.montecarlo_nu.value,
            last_line_interaction_angstrom=self.radial1d_model.last_line_interaction_angstrom.value,
            j_blues_norm_factor=self.radial1d_model.j_blues_norm_factor.value)


if __name__ == '__main__':
    # Quickly do the task the script it supposed to do if executed directly from terminal
    collector = Collector('abn_tom_test', 'abn_tom_test.yml',
                          'abn_tom_test_abundances.dat', 'abn_tom_test_densities.dat')
    collector.perform_run()
    collector.collect_astropy_quantities_from_radial1d_model()
    collector.collect_ndarrays_from_radial1d_model()
