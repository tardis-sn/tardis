#Module to read the rather complex config data
#Currently the configuration file is documented in tardis/data/example_configuration.ini

from astropy import constants, units
from ConfigParser import ConfigParser
import logging
import numpy as np
import os
import h5py
import pandas as pd

logger = logging.getLogger(__name__)

data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))

def calculate_w7_branch85_densities(velocities, time_explosion, time_0=19.9999584, density_coefficient=3e29):
    """
        Generated densities from the fit to W7 in Branch 85 page 620 (citation missing)

        Parameters
        ----------

        velocities : `~numpy.ndarray`
            velocities in cm/s

        time_explosion : `float`
            time since explosion needed to descale density with expansion

        time_0 : `float`
            time in seconds of the w7 model - default 19.999, no reason to change

        density_coefficient : `float`
            coefficient for the polynomial - obtained by fitting to W7, no reason to change

    """
    densities = density_coefficient * (velocities * 1e-5) ** -7
    densities *= (time_explosion / time_0) ** -3

    return densities[1:]


def read_w7_densities(fname=None):
    """
        Reading the density set for W7 in the density set h5 file

        Parameters
        ----------

        fname : `str`
            default None - defaults to tardis/data/density_sets.h5

    """
    pass


def read_lucy99_abundances(fname=None):
    """
    Reading the density set for W7 in the density set h5 file

    Parameters
    ----------

    fname : `str`
        default None - defaults to tardis/data/abundance_sets.h5
"""
    if fname is None:
        fname = os.path.join(data_dir, 'abundance_sets.h5')

    lucy99 = h5py.File(fname)['lucy99']

    logger.info("Choosing uniform abundance set 'lucy99':\n %s", pd.DataFrame(lucy99.__array__()))

    return dict(zip(lucy99.dtype.names, lucy99[0]))


def read_config(fname):
    config_object = ConfigParser()
    config_object.read(fname)
    tardis_configuration = TardisConfiguration()
    general_dict = dict(config_object.items('general'))
    parse_general_section(general_dict, tardis_configuration)
    abundance_dict = dict(config_object.items('abundances'))
    tardis_configuration.abundances = parse_abundance_section(abundance_dict)
    return tardis_configuration


class TardisConfiguration(object):
    pass


def parse_abundance_section(abundance_dict):
    abundance_set = abundance_dict.get('abundance_set', None)

    if abundance_set == 'lucy99':
        abundances = read_lucy99_abundances()
    else:
        raise ValueError('Currently only abundance_set=lucy99 supported')

    return abundances


def parse_general_section(config_dict, general_config):
    model_type = config_dict.pop('model_type')

    if model_type != 'radial1d':
        raise ValueError("Only supporting 'radial1d' at the moment")

    #reading time since explosion
    time_explosion_value, time_explosion_unit = config_dict.pop('time_explosion').split()
    general_config.time_explosion = units.Quantity(float(time_explosion_value), time_explosion_unit).to('s').value

    #Reading luminosity, special unit log_l_sun is luminosity given in log10 of solar units
    luminosity_value, luminosity_unit = config_dict.pop('luminosity').split()
    if luminosity_unit == 'log_lsun':
        general_config.luminosity = 10 ** (float(luminosity_value) + np.log10(constants.cgs.L_sun.value))
    else:
        general_config.luminosity = units.Quantity(float(luminosity_value), luminosity_unit)

    #reading number of shells
    no_of_shells = int(config_dict.pop('zones'))


    #reading velocities
    #set of velocities currently supported are v_inner, v_outer and v_sampling linear

    v_inner_value, v_inner_unit = config_dict.pop('v_inner').split()
    v_inner = units.Quantity(float(v_inner_value), v_inner_unit).to('cm/s').value

    v_outer_value, v_outer_unit = config_dict.pop('v_outer').split()
    v_outer = units.Quantity(float(v_outer_value), v_outer_unit).to('cm/s').value

    v_sampling = config_dict.pop('v_sampling')
    if v_sampling == 'linear':
        general_config.velocities = np.linspace(v_inner, v_outer, no_of_shells)
    else:
        raise ValueError('Currently only v_sampling = linear is possible')

    density_set = config_dict.pop('density_set')

    if density_set == 'w7_branch85':
        general_config.densities = calculate_w7_branch85_densities(general_config.velocities,
            general_config.time_explosion)
    else:
        raise ValueError('Curently only density_set = w7_branch85 is supported')


    #reading plasma type
    general_config.plasma_type = config_dict.pop('plasma_type')

    #reading line interaction type
    general_config.line_interaction_type = config_dict.pop('line_interaction_type')

    #reading number of packets and iterations
    general_config.calibration_packets = int(float(config_dict.pop('calibration_packets')))
    general_config.spectrum_packets = int(float(config_dict.pop('spectrum_packets')))
    general_config.iterations = int(float(config_dict.pop('iterations')))


    #TODO fix quantity spectral in astropy

    spectrum_start_value, spectrum_end_unit = config_dict.pop('spectrum_start').split()
    general_config.spectrum_start = units.Quantity(float(spectrum_start_value), spectrum_end_unit).value
    spectrum_end_value, spectrum_end_unit = config_dict.pop('spectrum_end').split()
    general_config.spectrum_end = units.Quantity(float(spectrum_end_value), spectrum_end_unit).value

    if config_dict != {}:
        logger.warn('Not all config options parsed - ignored %s' % config_dict)

    return general_config