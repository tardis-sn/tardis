"""
This is a simple script to collect the benchmark data by performing a single run of
provided setup. This benchmark data will be kept for future assertions. This script
is only to be run by trusted version of TARDIS which is certain to provide correct
output radial 1D model from the run.
"""

import os
import yaml
import numpy as np
from tardis import run_tardis


def data_path(model_dir, fname):
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
    return os.path.join(os.path.join(os.path.dirname(os.path.realpath(__file__)), model_dir, fname))


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


def perform_run(model_dir, config, abundances=None, densities=None):
    """Performs a single run of TARDIS with the specified model.

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

    Returns:
        expected_model: (~tardis.model.Radial1DModel)
            The obtained model from the TARDIS run with specified setup.
    """
    # Path to the yaml config files and atom data files
    yml_config_filepath = data_path(model_dir, config)
    save_kurucz_atom_data_file()

    # Obtain a dictionary out of yaml config filepath and override values of the
    # path to atom data, densities and/or abundances data files if provided
    config_dict = yaml.load(open(yml_config_filepath))
    config_dict['atom_data'] = os.path.abspath('/tmp/kurucz_cd23_chianti_H_He.h5')
    if abundances is not None:
        config_dict['model']['abundances']['filename'] = data_path(model_dir, abundances)

    if densities is not None:
        config_dict['model']['structure']['filename'] = data_path(model_dir, densities)

    # Prepare a standard model from a 'believed' stable build of TARDIS
    expected_model = run_tardis(config_dict)
    return expected_model


def collect_ndarrays_from_radial1d_model(expected_model):
    # the obtained 'expected_model' contains certain members as numpy arrays
    # all are collected and put into one compressed 'expected_ndarrays.npz' file
    # this 'expected_ndarrays.npz' can be used for subsequent assertions in testing
    np.savez_compressed(
            data_path('expected_ndarrays.npz'),
            last_interaction_type=expected_model.last_interaction_type,
            last_line_interaction_out_id=expected_model.last_line_interaction_out_id,
            last_line_interaction_in_id=expected_model.last_line_interaction_in_id,
            j_estimators=expected_model.j_estimators,
            j_blue_estimators=expected_model.j_blue_estimators,
            last_line_interaction_shell_id=expected_model.last_line_interaction_shell_id,
            nubar_estimators=expected_model.nubar_estimators,
            ws=expected_model.ws)


def collect_astropy_quantities_from_radial1d_model(expected_model):
    # apart from the ndarrays, there are certain 'astropy.units.quantity.Quantity' type
    # members. their values are numpy ndarrays so are collected and put together in a file
    # 'expected_quantities.npz' - this can be used for subsequent assertions in testing
    np.savez_compressed(
            data_path('expected_quantities.npz'),
            t_rads=expected_model.t_rads.value,
            luminosity_inner=expected_model.luminosity_inner.value,
            montecarlo_luminosity=expected_model.montecarlo_luminosity.value,
            montecarlo_virtual_luminosity=expected_model.montecarlo_virtual_luminosity.value,
            time_of_simulation=expected_model.time_of_simulation.value,
            montecarlo_nu=expected_model.montecarlo_nu.value,
            last_line_interaction_angstrom=expected_model.last_line_interaction_angstrom.value,
            j_blues_norm_factor=expected_model.j_blues_norm_factor.value)


if __name__ == '__main__':
    # Quickly do the task the script it supposed to do if executed directly from terminal
    mdl = collect_expected_data()
