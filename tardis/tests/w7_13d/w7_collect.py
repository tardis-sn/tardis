"""
This is a simple script to collect the data from a single run of Stratified W7 model.
The data is meant to be set as a benchmark for future assertions. This script is only
to be run with a trusted version of TARDIS which is certain to provide correct output
radial 1D model from the run.
"""

import os
import numpy as np
from tardis import run_tardis


def data_path(fname):
    # $TARDIS_PATH/tests/data/w7_13d/'fname'
    return os.path.join(os.path.join(os.path.dirname(os.path.realpath(__file__)), fname))


def collect_expected_data():
    # path to the yaml config files and atom data files
    yml_config_file = data_path('tardis_w7.yml')

    # to keep the tardis working repo clean from a large atom data file, we store
    # this file in the 'tmp' directory, and notify the user about its existence
    if not os.path.exists(os.path.abspath('/tmp/kurucz_cd23_chianti_H_He.h5')):
        os.system('wget http://www.mpa-garching.mpg.de/~michi/tardis/data/kurucz_cd23_chianti_H_He.zip')
        os.system('unzip kurucz_cd23_chianti_H_He.zip')
        os.system('mv kurucz_cd23_chianti_H_He.h5 /tmp/')
        os.system('rm kurucz_cd23_chianti_H_He.zip')
        print('Successfully obtained atom data file and store in /tmp directory.')
    else:
        print('The atom data file already exists in /tmp directory.')

    atom_data_file = os.path.abspath('/tmp/kurucz_cd23_chianti_H_He.h5')

    # prepare a standard model from a 'believed' stable build of TARDIS
    # this is a long run and takes upto 30 to 35 minutes
    expected_model = run_tardis(config=yml_config_file, atom_data=atom_data_file)

    # collect properties stored as ndarrays, and 'astropy.units.quantity.Quantity's
    # these would be stored in two separate files in compressed binary formats and
    # they will be used as a benchmark while doing a run for testing
    collect_ndarrays_from_radial1d_model(expected_model)
    collect_astropy_quantities_from_radial1d_model(expected_model)

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
