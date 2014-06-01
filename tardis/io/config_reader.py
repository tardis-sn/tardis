# Module to read the rather complex config data

import logging
import os
import pprint

from astropy import constants, units as u
import numpy as np
import pandas as pd
import yaml

import tardis
from tardis.io.model_reader import read_density_file, \
    calculate_density_after_time, read_abundances_file
from tardis.io.config_validator import Config
from tardis import atomic
from tardis.util import species_string_to_tuple, parse_quantity, \
    element_symbol2atomic_number

import copy

pp = pprint.PrettyPrinter(indent=4)

logger = logging.getLogger(__name__)

data_dir = os.path.join(tardis.__path__[0], 'data')

default_config_definition_file = os.path.join(data_dir,
                                              'tardis_config_definition.yml')
#File parsers for different file formats:


density_structure_fileparser = {}

inv_ni56_efolding_time = 1 / (8.8 * u.day)
inv_co56_efolding_time = 1 / (113.7 * u.day)
inv_cr48_efolding_time = 1 / (1.29602 * u.day)
inv_v48_efolding_time = 1 / (23.0442 * u.day)
inv_fe52_efolding_time = 1 / (0.497429 * u.day)
inv_mn52_efolding_time = 1 / (0.0211395 * u.day)


class ConfigurationError(ValueError):
    pass


def parse_quantity_linspace(quantity_linspace_dictionary, add_one=True):
    """
    parse a dictionary of the following kind
    {'start': 5000 km/s,
     'stop': 10000 km/s,
     'num': 1000}

    Parameters
    ----------

    quantity_linspace_dictionary: ~dict

    add_one: boolean, default: True

    Returns
    -------

    ~np.array

    """

    start = parse_quantity(quantity_linspace_dictionary['start'])
    stop = parse_quantity(quantity_linspace_dictionary['stop'])

    try:
        stop = stop.to(start.unit)
    except u.UnitsError:
        raise ConfigurationError('"start" and "stop" keyword must be compatible quantities')

    num = quantity_linspace_dictionary['num']
    if add_one:
        num += 1

    return np.linspace(start.value, stop.value, num=num) * start.unit


def parse_spectral_bin(spectral_bin_boundary_1, spectral_bin_boundary_2):
    spectral_bin_boundary_1 = parse_quantity(spectral_bin_boundary_1).to('Angstrom', u.spectral())
    spectral_bin_boundary_2 = parse_quantity(spectral_bin_boundary_2).to('Angstrom', u.spectral())

    spectrum_start_wavelength = min(spectral_bin_boundary_1, spectral_bin_boundary_2)
    spectrum_end_wavelength = max(spectral_bin_boundary_1, spectral_bin_boundary_2)

    return spectrum_start_wavelength, spectrum_end_wavelength


def calculate_exponential_density(velocities, v_0, rho0):
    """
    This function computes the exponential density profile.
    :math:`\\rho = \\rho_0 \\times \\exp \\left( -\\frac{v}{v_0} \\right)`

    Parameters
    ----------

    velocities : ~astropy.Quantity
        Array like velocity profile
    velocity_0 : ~astropy.Quantity
        reference velocity
    rho0 : ~astropy.Quantity
        reference density

    Returns
    -------

    densities : ~astropy.Quantity

    """
    densities = rho0 * np.exp(-(velocities / v_0))
    return densities


def calculate_power_law_density(velocities, velocity_0, rho_0, exponent):
    """

    This function computes a descret exponential density profile.
    :math:`\\rho = \\rho_0 \\times \\left( \\frac{v}{v_0} \\right)^n`

    Parameters
    ----------

    velocities : ~astropy.Quantity
        Array like velocity profile
    velocity_0 : ~astropy.Quantity
        reference velocity
    rho0 : ~astropy.Quantity
        reference density
    exponent : ~float
        exponent used in the powerlaw

    Returns
    -------

    densities : ~astropy.Quantity

    """
    densities = rho_0 * np.power((velocities / velocity_0), exponent)
    return densities


def parse_model_file_section(model_setup_file_dict, time_explosion):
    def parse_artis_model_setup_files(model_file_section_dict, time_explosion):

        ###### Reading the structure part of the ARTIS file pair
        structure_fname = model_file_section_dict['structure_fname']

        for i, line in enumerate(file(structure_fname)):
            if i == 0:
                no_of_shells = np.int64(line.strip())
            elif i == 1:
                time_of_model = u.Quantity(float(line.strip()), 'day').to('s')
            elif i == 2:
                break

        artis_model_columns = ['velocities', 'mean_densities_0', 'ni56_fraction', 'co56_fraction', 'fe52_fraction',
                               'cr48_fraction']
        artis_model = np.recfromtxt(structure_fname, skip_header=2, usecols=(1, 2, 4, 5, 6, 7), unpack=True,
                                    dtype=[(item, np.float64) for item in artis_model_columns])
        #converting densities from log(g/cm^3) to g/cm^3 and stretching it to the current ti
        velocities = u.Quantity(np.append([0], artis_model['velocities']), 'km/s').to('cm/s')
        mean_densities_0 = u.Quantity(10 ** artis_model['mean_densities_0'], 'g/cm^3')

        mean_densities = calculate_density_after_time(mean_densities_0, time_of_model, time_explosion)


        #Verifying information
        if len(mean_densities) == no_of_shells:
            logger.debug('Verified ARTIS model structure file %s (no_of_shells=length of dataset)', structure_fname)
        else:
            raise ConfigurationError(
                'Error in ARTIS file %s - Number of shells not the same as dataset length' % structure_fname)

        v_inner = velocities[:-1]
        v_outer = velocities[1:]

        volumes = (4 * np.pi / 3) * (time_of_model ** 3) * ( v_outer ** 3 - v_inner ** 3)
        masses = (volumes * mean_densities_0 / constants.M_sun).to(1)

        logger.info('Read ARTIS configuration file %s - found %d zones with total mass %g Msun', structure_fname,
                    no_of_shells, sum(masses.value))

        if 'v_lowest' in model_file_section_dict:
            v_lowest = parse_quantity(model_file_section_dict['v_lowest']).to('cm/s').value
            min_shell = v_inner.value.searchsorted(v_lowest)
        else:
            min_shell = 1

        if 'v_highest' in model_file_section_dict:
            v_highest = parse_quantity(model_file_section_dict['v_highest']).to('cm/s').value
            max_shell = v_outer.value.searchsorted(v_highest)
        else:
            max_shell = no_of_shells
        artis_model = artis_model[min_shell:max_shell]
        v_inner = v_inner[min_shell:max_shell]
        v_outer = v_outer[min_shell:max_shell]
        mean_densities = mean_densities[min_shell:max_shell]

        ###### Reading the abundance part of the ARTIS file pair
        abundances_fname = model_file_section_dict['abundances_fname']
        abundances = pd.DataFrame(np.loadtxt(abundances_fname)[min_shell:max_shell, 1:].transpose(),
                                  index=np.arange(1, 31))

        ni_stable = abundances.ix[28] - artis_model['ni56_fraction']
        co_stable = abundances.ix[27] - artis_model['co56_fraction']
        fe_stable = abundances.ix[26] - artis_model['fe52_fraction']
        mn_stable = abundances.ix[25] - 0.0
        cr_stable = abundances.ix[24] - artis_model['cr48_fraction']
        v_stable = abundances.ix[23] - 0.0
        ti_stable = abundances.ix[22] - 0.0

        abundances.ix[28] = ni_stable
        abundances.ix[28] += artis_model['ni56_fraction'] * np.exp(
            -(time_explosion * inv_ni56_efolding_time).to(1).value)

        abundances.ix[27] = co_stable
        abundances.ix[27] += artis_model['co56_fraction'] * np.exp(
            -(time_explosion * inv_co56_efolding_time).to(1).value)
        abundances.ix[27] += (inv_ni56_efolding_time * artis_model['ni56_fraction'] /
                              (inv_ni56_efolding_time - inv_co56_efolding_time)) * \
                             (np.exp(-(inv_co56_efolding_time * time_explosion).to(1).value) - np.exp(
                                 -(inv_ni56_efolding_time * time_explosion).to(1).value))

        abundances.ix[26] = fe_stable
        abundances.ix[26] += artis_model['fe52_fraction'] * np.exp(
            -(time_explosion * inv_fe52_efolding_time).to(1).value)
        abundances.ix[26] += ((artis_model['co56_fraction'] * inv_ni56_efolding_time
                               - artis_model['co56_fraction'] * inv_co56_efolding_time
                               + artis_model['ni56_fraction'] * inv_ni56_efolding_time
                               - artis_model['ni56_fraction'] * inv_co56_efolding_time
                               - artis_model['co56_fraction'] * inv_ni56_efolding_time * np.exp(
            -(inv_co56_efolding_time * time_explosion).to(1).value)
                               + artis_model['co56_fraction'] * inv_co56_efolding_time * np.exp(
            -(inv_co56_efolding_time * time_explosion).to(1).value)
                               - artis_model['ni56_fraction'] * inv_ni56_efolding_time * np.exp(
            -(inv_co56_efolding_time * time_explosion).to(1).value)
                               + artis_model['ni56_fraction'] * inv_co56_efolding_time * np.exp(
            -(inv_ni56_efolding_time * time_explosion).to(1).value))
                              / (inv_ni56_efolding_time - inv_co56_efolding_time))

        abundances.ix[25] = mn_stable
        abundances.ix[25] += (inv_fe52_efolding_time * artis_model['fe52_fraction'] /
                              (inv_fe52_efolding_time - inv_mn52_efolding_time)) * \
                             (np.exp(-(inv_mn52_efolding_time * time_explosion).to(1).value) - np.exp(
                                 -(inv_fe52_efolding_time * time_explosion).to(1).value))

        abundances.ix[24] = cr_stable
        abundances.ix[24] += artis_model['cr48_fraction'] * np.exp(
            -(time_explosion * inv_cr48_efolding_time).to(1).value)
        abundances.ix[24] += ((artis_model['fe52_fraction'] * inv_fe52_efolding_time
                               - artis_model['fe52_fraction'] * inv_mn52_efolding_time
                               - artis_model['fe52_fraction'] * inv_fe52_efolding_time * np.exp(
            -(inv_mn52_efolding_time * time_explosion).to(1).value)
                               + artis_model['fe52_fraction'] * inv_mn52_efolding_time * np.exp(
            -(inv_fe52_efolding_time * time_explosion).to(1).value))
                              / (inv_fe52_efolding_time - inv_mn52_efolding_time))

        abundances.ix[23] = v_stable
        abundances.ix[23] += (inv_cr48_efolding_time * artis_model['cr48_fraction'] /
                              (inv_cr48_efolding_time - inv_v48_efolding_time)) * \
                             (np.exp(-(inv_v48_efolding_time * time_explosion).to(1).value) - np.exp(
                                 -(inv_cr48_efolding_time * time_explosion).to(1).value))

        abundances.ix[22] = ti_stable
        abundances.ix[22] += ((artis_model['cr48_fraction'] * inv_cr48_efolding_time
                               - artis_model['cr48_fraction'] * inv_v48_efolding_time
                               - artis_model['cr48_fraction'] * inv_cr48_efolding_time * np.exp(
            -(inv_v48_efolding_time * time_explosion).to(1).value)
                               + artis_model['cr48_fraction'] * inv_v48_efolding_time * np.exp(
            -(inv_cr48_efolding_time * time_explosion).to(1).value))
                              / (inv_cr48_efolding_time - inv_v48_efolding_time))

        if 'split_shells' in model_file_section_dict:
            split_shells = int(model_file_section_dict['split_shells'])
        else:
            split_shells = 1

        if split_shells > 1:
            logger.info('Increasing the number of shells by a factor of %s' % split_shells)
            no_of_shells = len(v_inner)
            velocities = np.linspace(v_inner[0], v_outer[-1], no_of_shells * split_shells + 1)
            v_inner = velocities[:-1]
            v_outer = velocities[1:]
            old_mean_densities = mean_densities
            mean_densities = np.empty(no_of_shells * split_shells) * old_mean_densities.unit
            new_abundance_data = np.empty((abundances.values.shape[0], no_of_shells * split_shells))
            for i in xrange(split_shells):
                mean_densities[i::split_shells] = old_mean_densities
                new_abundance_data[:, i::split_shells] = abundances.values

            abundances = pd.DataFrame(new_abundance_data, index=abundances.index)




            #def parser_simple_ascii_model

        return v_inner, v_outer, mean_densities, abundances

    model_file_section_parser = {}
    model_file_section_parser['artis'] = parse_artis_model_setup_files

    try:
        parser = model_file_section_parser[model_setup_file_dict['type']]
    except KeyError:
        raise ConfigurationError('In abundance file section only types %s are allowed (supplied %s) ' %
                                 (model_file_section_parser.keys(), model_file_section_parser['type']))

    return parser(model_setup_file_dict, time_explosion)


def parse_density_file_section(density_file_dict, time_explosion):
    density_file_parser = {}

    def parse_artis_density(density_file_dict, time_explosion):
        density_file = density_file_dict['name']
        for i, line in enumerate(file(density_file)):
            if i == 0:
                no_of_shells = np.int64(line.strip())
            elif i == 1:
                time_of_model = u.Quantity(float(line.strip()), 'day').to('s')
            elif i == 2:
                break

        velocities, mean_densities_0 = np.recfromtxt(density_file, skip_header=2, usecols=(1, 2), unpack=True)
        #converting densities from log(g/cm^3) to g/cm^3 and stretching it to the current ti
        velocities = u.Quantity(np.append([0], velocities), 'km/s').to('cm/s')
        mean_densities_0 = u.Quantity(10 ** mean_densities_0, 'g/cm^3')

        mean_densities = calculate_density_after_time(mean_densities_0, time_of_model, time_explosion)


        #Verifying information
        if len(mean_densities) == no_of_shells:
            logger.debug('Verified ARTIS file %s (no_of_shells=length of dataset)', density_file)
        else:
            raise ConfigurationError(
                'Error in ARTIS file %s - Number of shells not the same as dataset length' % density_file)

        min_shell = 1
        max_shell = no_of_shells

        v_inner = velocities[:-1]
        v_outer = velocities[1:]

        volumes = (4 * np.pi / 3) * (time_of_model ** 3) * ( v_outer ** 3 - v_inner ** 3)
        masses = (volumes * mean_densities_0 / constants.M_sun).to(1)

        logger.info('Read ARTIS configuration file %s - found %d zones with total mass %g Msun', density_file,
                    no_of_shells, sum(masses.value))

        if 'v_lowest' in density_file_dict:
            v_lowest = parse_quantity(density_file_dict['v_lowest']).to('cm/s').value
            min_shell = v_inner.value.searchsorted(v_lowest)
        else:
            min_shell = 1

        if 'v_highest' in density_file_dict:
            v_highest = parse_quantity(density_file_dict['v_highest']).to('cm/s').value
            max_shell = v_outer.value.searchsorted(v_highest)
        else:
            max_shell = no_of_shells

        v_inner = v_inner[min_shell:max_shell]
        v_outer = v_outer[min_shell:max_shell]
        mean_densities = mean_densities[min_shell:max_shell]

        return v_inner, v_outer, mean_densities, min_shell, max_shell

    density_file_parser['artis'] = parse_artis_density

    try:
        parser = density_file_parser[density_file_dict['type']]
    except KeyError:
        raise ConfigurationError('In abundance file section only types %s are allowed (supplied %s) ' %
                                 (density_file_parser.keys(), density_file_dict['type']))

    return parser(density_file_dict, time_explosion)


def parse_density_section(density_dict, v_inner, v_outer, time_explosion):
    density_parser = {}


    #Parse density uniform
    def parse_uniform(density_dict, v_inner, v_outer, time_explosion):
        no_of_shells = len(v_inner)
        return parse_quantity(density_dict['value']).to('g cm^-3') * \
               np.ones(no_of_shells)

    density_parser['uniform'] = parse_uniform

    #Parse density branch85 w7
    def parse_branch85(density_dict, v_inner, v_outer, time_explosion):
        velocities = 0.5 * (v_inner + v_outer)

        densities = calculate_power_law_density(velocities,
                                                density_dict['w7_v_0'],
                                                density_dict['w7_rho_0'], -7)

        densities = calculate_density_after_time(densities,
                                                 density_dict['w7_time_0'],
                                                 time_explosion)


        return densities

    density_parser['branch85_w7'] = parse_branch85

    def parse_power_law(density_dict, v_inner, v_outer, time_explosion):
        time_0 = density_dict.pop('time_0')
        rho_0 = density_dict.pop('rho_0')
        v_0 = density_dict.pop('v_0')
        exponent = density_dict.pop('exponent')
        velocities = 0.5 * (v_inner + v_outer)
        densities = calculate_power_law_density(velocities, v_0, rho_0, exponent)
        densities = calculate_density_after_time(densities, time_0, time_explosion)
        return densities

    density_parser['power_law'] = parse_power_law

    def parse_exponential(density_dict, v_inner, v_outer, time_explosion):
        time_0 = density_dict.pop('time_0')
        rho_0 = density_dict.pop('rho_0')
        v_0 = density_dict.pop('v_0')

        velocities = 0.5 * (v_inner + v_outer)
        densities = calculate_exponential_density(velocities, v_0, rho_0)
        densities = calculate_density_after_time(densities, time_0, time_explosion)
        return densities

    density_parser['exponential'] = parse_exponential

    try:
        parser = density_parser[density_dict['type']]
    except KeyError:
        raise ConfigurationError('In density section only types %s are allowed (supplied %s) ' %
                                 (density_parser.keys(), density_dict['type']))
    return parser(density_dict, v_inner, v_outer, time_explosion)


def parse_abundance_file_section(abundance_file_dict, abundances, min_shell, max_shell):
    abundance_file_parser = {}

    def parse_artis(abundance_file_dict, abundances, min_shell, max_shell):
        #### ---- debug ----
        time_of_model = 0.0

        ####
        fname = abundance_file_dict['name']
        max_atom = 30
        logger.info("Parsing ARTIS Abundance section from shell %d to %d", min_shell, max_shell)

        abundances.values[:max_atom, :] = np.loadtxt(fname)[min_shell:max_shell, 1:].transpose()

        return abundances

    abundance_file_parser['artis'] = parse_artis

    try:
        parser = abundance_file_parser[abundance_file_dict['type']]
    except KeyError:
        raise ConfigurationError('In abundance file section only types %s are allowed (supplied %s) ' %
                                 (abundance_file_parser.keys(), abundance_file_dict['type']))

    return parser(abundance_file_dict, abundances, min_shell, max_shell)


def parse_supernova_section(supernova_dict):
    """
    Parse the supernova section

    Parameters
    ----------

    supernova_dict: dict
        YAML parsed supernova dict

    Returns
    -------

    config_dict: dict

    """
    config_dict = {}

    #parse luminosity
    luminosity_value, luminosity_unit = supernova_dict['luminosity_requested'].strip().split()

    if luminosity_unit == 'log_lsun':
        config_dict['luminosity_requested'] = 10 ** (
        float(luminosity_value) + np.log10(constants.L_sun.cgs.value)) * u.erg / u.s
    else:
        config_dict['luminosity_requested'] = (float(luminosity_value) * u.Unit(luminosity_unit)).to('erg/s')

    config_dict['time_explosion'] = parse_quantity(supernova_dict['time_explosion']).to('s')

    if 'distance' in supernova_dict:
        config_dict['distance'] = parse_quantity(supernova_dict['distance'])
    else:
        config_dict['distance'] = None

    if 'luminosity_wavelength_start' in supernova_dict:
        config_dict['luminosity_nu_end'] = parse_quantity(supernova_dict['luminosity_wavelength_start']). \
            to('Hz', u.spectral())
    else:
        config_dict['luminosity_nu_end'] = np.inf * u.Hz

    if 'luminosity_wavelength_end' in supernova_dict:
        config_dict['luminosity_nu_start'] = parse_quantity(supernova_dict['luminosity_wavelength_end']). \
            to('Hz', u.spectral())
    else:
        config_dict['luminosity_nu_start'] = 0.0 * u.Hz

    return config_dict


def parse_spectrum_list2dict(spectrum_list):
    """
    Parse the spectrum list [start, stop, num] to a list
    """

    if spectrum_list[0].unit.physical_type != 'length' and \
                    spectrum_list[1].unit.physical_type != 'length':
        raise ValueError('start and end of spectrum need to be a length')


    spectrum_config_dict = {}
    spectrum_config_dict['start'] = spectrum_list[0]
    spectrum_config_dict['end'] = spectrum_list[1]
    spectrum_config_dict['bins'] = spectrum_list[2]

    spectrum_frequency = np.linspace(
        spectrum_config_dict['end'].to('Hz', u.spectral()),
        spectrum_config_dict['start'].to('Hz', u.spectral()),
        num=spectrum_config_dict['bins'] + 1)

    spectrum_config_dict['frequency'] = spectrum_frequency

    return spectrum_config_dict



def parse_convergence_section(convergence_section_dict):
    """
    Parse the convergence section dictionary

    Parameters
    ----------

    convergence_section_dict: ~dict
        dictionary
    """



    for convergence_variable in ['t_inner', 't_rad', 'w']:
        if convergence_variable not in convergence_section_dict:
            convergence_section_dict[convergence_variable] = {}

        updated_convergence_dict = convergence_section_dict[
            'global_convergence_parameters'].copy()
        updated_convergence_dict.update(
            convergence_section_dict[convergence_variable])

        convergence_section_dict[convergence_variable] = \
            updated_convergence_dict

    return convergence_section_dict

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
    densities = calculate_density_after_time(densities, time_0, time_explosion)

    return densities[1:]


class ConfigurationNameSpace(dict):
    """
    The configuration name space class allows to wrap a dictionary and adds
    utility functions for easy access. Accesses like a.b.c are then possible

    Code from http://goo.gl/KIaq8I

    Parameters
    ----------

    config_dict: ~dict
        configuration dictionary

    Returns
    -------

    config_ns: ConfigurationNameSpace
    """

    marker = object()
    def __init__(self, value=None):
        if value is None:
            pass
        elif isinstance(value, dict):
            for key in value:
                self.__setitem__(key, value[key])
        else:
            raise TypeError, 'expected dict'

    def __setitem__(self, key, value):
        if isinstance(value, dict) and not isinstance(value,
                                                      ConfigurationNameSpace):
            value = ConfigurationNameSpace(value)

        if key in self and hasattr(self[key], 'unit'):
            value = u.Quantity(value, self[key].unit)

        dict.__setitem__(self, key, value)

    def __getitem__(self, key):
        return super(ConfigurationNameSpace, self).__getitem__(key)

    def __getattr__(self, item):
        if item in self:
            return self[item]
        else:
            super(ConfigurationNameSpace, self).__getattribute__(item)

    __setattr__ = __setitem__

    def __dir__(self):
        self.keys()

    def get_config_item(self, config_item_string):
        """
        Get configuration items using a string of type 'a.b.param'

        Parameters
        ----------

        config_item_string: ~str
            string of shape 'section1.sectionb.param1'
        """
        config_item_path = config_item_string.split('.')

        if len(config_item_path) == 1:
            config_item = config_item_path[0]

            if config_item.startswith('item'):
                return self[config_item_path[0]]
            else:
                return self[config_item]
        elif len(config_item_path) == 2 and\
                config_item_path[1].startswith('item'):
            return self[config_item_path[0]][
                int(config_item_path[1].replace('item', ''))]

        else:
            return self[config_item_path[0]].get_config_item(
                '.'.join(config_item_path[1:]))

    def set_config_item(self, config_item_string, value):
        """
        set configuration items using a string of type 'a.b.param'

        Parameters
        ----------

        config_item_string: ~str
            string of shape 'section1.sectionb.param1'

        value:
            value to set the parameter with it
        """

        config_item_path = config_item_string.split('.')
        if len(config_item_path) == 1:
            self[config_item_path[0]] = value
        elif len(config_item_path) == 2 and \
                config_item_path[1].startswith('item'):
            current_value = self[config_item_path[0]][
                int(config_item_path[1].replace('item', ''))]
            if hasattr(current_value, 'unit'):
                self[config_item_path[0]][
                    int(config_item_path[1].replace('item', ''))] =\
                    u.Quantity(value, current_value.unit)
            else:
                self[config_item_path[0]][
                    int(config_item_path[1].replace('item', ''))] = value

        else:
            self[config_item_path[0]].set_config_item(
                '.'.join(config_item_path[1:]), value)

    def deepcopy(self):
        return ConfigurationNameSpace(copy.deepcopy(dict(self)))


class Configuration(ConfigurationNameSpace):
    """
    Tardis configuration class
    """

    @classmethod
    def from_yaml(cls, fname, test_parser=False):
        try:
            yaml_dict = yaml.load(file(fname))
        except IOError as e:
            logger.critical('No config file named: %s', fname)
            raise e

        tardis_config_version = yaml_dict.get('tardis_config_version', None)
        if tardis_config_version != 'v1.0':
            raise ConfigurationError('Currently only tardis_config_version v1.0 supported')

        return cls.from_config_dict(yaml_dict, test_parser=test_parser)

    @classmethod
    def from_config_dict(cls, config_dict, atom_data=None, test_parser=False,
                         config_definition_file=None):
        """
        Validating and subsequently parsing a config file.


        Parameters
        ----------

        config_dict : ~dict
            dictionary of a raw unvalidated config file

        atom_data: ~tardis.atomic.AtomData
            atom data object. if `None` will be tried to be read from
            atom data file path in the config_dict [default=None]

        test_parser: ~bool
            switch on to ignore a working atom_data, mainly useful for
            testing this reader

        config_definition_file: ~str
            path to config definition file, if `None` will be set to the default
            in the `data` directory that ships with TARDIS


        Returns
        -------

        `tardis.config_reader.Configuration`

        """

        if config_definition_file is None:
            config_definition_file = default_config_definition_file

        config_definition = yaml.load(open(config_definition_file))

        validated_config_dict = Config(config_definition,
                                       config_dict).get_config()

        #First let's see if we can find an atom_db anywhere:
        if test_parser:
            atom_data = None
        elif 'atom_data' in validated_config_dict.keys():
            atom_data_fname = validated_config_dict['atom_data']
            validated_config_dict['atom_data_fname'] = atom_data_fname
        else:
            raise ConfigurationError('No atom_data key found in config or command line')

        if atom_data is None and not test_parser:
            logger.info('Reading Atomic Data from %s', atom_data_fname)
            atom_data = atomic.AtomData.from_hdf5(atom_data_fname)
        else:
            atom_data = atom_data



        #Parsing supernova dictionary
        validated_config_dict['supernova']['luminosity_nu_start'] = \
            validated_config_dict['supernova']['luminosity_wavelength_end'].to(
                u.Hz, u.spectral())
        try:
            validated_config_dict['supernova']['luminosity_nu_end'] = \
                (validated_config_dict['supernova']
                 ['luminosity_wavelength_start'].to(u.Hz, u.spectral()))
        except ZeroDivisionError:
            validated_config_dict['supernova']['luminosity_nu_end'] = (
                np.inf * u.Hz)

        validated_config_dict['supernova']['time_explosion'] = (
            validated_config_dict['supernova']['time_explosion'].cgs)

        validated_config_dict['supernova']['luminosity_requested'] = (
            validated_config_dict['supernova']['luminosity_requested'].cgs)

        #Parsing the model section

        model_section = validated_config_dict['model']
        v_inner = None
        v_outer = None
        mean_densities = None
        abundances = None



        structure_section = model_section['structure']

        if structure_section['type'] == 'specific':
            start, stop, num = model_section['structure']['velocity']
            num += 1
            velocities = np.linspace(start, stop, num)

            v_inner, v_outer = velocities[:-1], velocities[1:]
            mean_densities = parse_density_section(
                model_section['structure']['density'], v_inner, v_outer,
                validated_config_dict['supernova']['time_explosion']).cgs

        elif structure_section['type'] == 'file':
            v_inner, v_outer, mean_densities, inner_boundary_index, \
            outer_boundary_index = read_density_file(
                structure_section['filename'], structure_section['filetype'],
                validated_config_dict['supernova']['time_explosion'],
                structure_section['v_inner_boundary'],
                structure_section['v_outer_boundary'])

        r_inner = validated_config_dict['supernova']['time_explosion'] * v_inner
        r_outer = validated_config_dict['supernova']['time_explosion'] * v_outer
        r_middle = 0.5 * (r_inner + r_outer)

        structure_validated_config_dict = {}
        structure_section['v_inner'] = v_inner.cgs
        structure_section['v_outer'] = v_outer.cgs
        structure_section['mean_densities'] = mean_densities.cgs
        no_of_shells = len(v_inner)
        structure_section['no_of_shells'] = no_of_shells
        structure_section['r_inner'] = r_inner.cgs
        structure_section['r_outer'] = r_outer.cgs
        structure_section['r_middle'] = r_middle.cgs
        structure_section['volumes'] = ((4. / 3) * np.pi * \
                                       (r_outer ** 3 -
                                        r_inner ** 3)).cgs


        #### TODO the following is legacy code and should be removed
        validated_config_dict['structure'] = \
            validated_config_dict['model']['structure']
        # ^^^^^^^^^^^^^^^^


        abundances_section = model_section['abundances']

        if abundances_section['type'] == 'uniform':
            abundances = pd.DataFrame(columns=np.arange(no_of_shells),
                                      index=pd.Index(np.arange(1, 120), name='atomic_number'), dtype=np.float64)

            for element_symbol_string in abundances_section:
                if element_symbol_string == 'type': continue
                z = element_symbol2atomic_number(element_symbol_string)
                abundances.ix[z] = float(abundances_section[element_symbol_string])

        elif abundances_section['type'] == 'file':
            index, abundances = read_abundances_file(abundances_section['filename'], abundances_section['filetype'],
                                                     inner_boundary_index, outer_boundary_index)
            if len(index) != no_of_shells:
                raise ConfigurationError('The abundance file specified has not the same number of cells'
                                         'as the specified density profile')

        abundances = abundances.replace(np.nan, 0.0)

        abundances = abundances[abundances.sum(axis=1) > 0]

        norm_factor = abundances.sum(axis=0)

        if np.any(np.abs(norm_factor - 1) > 1e-12):
            logger.warning("Abundances have not been normalized to 1. - normalizing")
            abundances /= norm_factor

        validated_config_dict['abundances'] = abundances



        ########### DOING PLASMA SECTION ###############
        plasma_section = validated_config_dict['plasma']

        if plasma_section['initial_t_inner'] < 0.0 * u.K:
            luminosity_requested = validated_config_dict['supernova']['luminosity_requested']
            plasma_section['t_inner'] = ((luminosity_requested /
                                          (4 * np.pi * r_inner[0] ** 2 *
                                           constants.sigma_sb)) ** .25).to('K')
            logger.info('"initial_t_inner" is not specified in the plasma '
                        'section - initializing to %s with given luminosity',
                        plasma_section['t_inner'])

        plasma_section['t_rads'] = np.ones(no_of_shells) * \
                                   plasma_section['initial_t_rad']
        if plasma_section['disable_electron_scattering'] is False:
            logger.debug("Electron scattering switched on")
            validated_config_dict['montecarlo']['sigma_thomson'] = 6.652486e-25 / (u.cm ** 2)
        else:
            logger.warn('Disabling electron scattering - this is not physical')
            validated_config_dict['montecarlo']['sigma_thomson'] = 1e-200 / (u.cm ** 2)




        ##### NLTE subsection of Plasma start
        nlte_validated_config_dict = {}
        nlte_species = []
        nlte_section = plasma_section['nlte']

        nlte_species_list = nlte_section.pop('species')
        for species_string in nlte_species_list:
            nlte_species.append(species_string_to_tuple(species_string))

        nlte_validated_config_dict['species'] = nlte_species
        nlte_validated_config_dict['species_string'] = nlte_species_list
        nlte_validated_config_dict.update(nlte_section)

        if 'coronal_approximation' not in nlte_section:
            logger.debug('NLTE "coronal_approximation" not specified in NLTE section - defaulting to False')
            nlte_validated_config_dict['coronal_approximation'] = False

        if 'classical_nebular' not in nlte_section:
            logger.debug('NLTE "classical_nebular" not specified in NLTE section - defaulting to False')
            nlte_validated_config_dict['classical_nebular'] = False


        elif nlte_section:  #checks that the dictionary is not empty
            logger.warn('No "species" given - ignoring other NLTE options given:\n%s',
                        pp.pformat(nlte_section))

        if not nlte_validated_config_dict:
            nlte_validated_config_dict['species'] = []

        plasma_section['nlte'] = nlte_validated_config_dict

        #^^^^^^^^^^^^^^ End of Plasma Section

        ##### Monte Carlo Section

        montecarlo_section = validated_config_dict['montecarlo']

        default_convergence_section = {'type': 'damped',
                                      'lock_t_inner_cyles': 1,
                                      't_inner_update_exponent': -0.5,
                                      'global_convergence_parameters' : {
                                          'damping_constant': 0.5}}



        if montecarlo_section['convergence_strategy'] is None:
            logger.warning('No convergence criteria selected - just damping by 0.5 for w, t_rad and t_inner')
            montecarlo_section['convergence_strategy'] = default_convergence_section
        else:
            1/0

        montecarlo_section['convergence_strategy'] = parse_convergence_section(
            montecarlo_section['convergence_strategy'])

        black_body_section = montecarlo_section['black_body_sampling']
        montecarlo_section['black_body_sampling'] = {}
        montecarlo_section['black_body_sampling']['start'] = \
            black_body_section[0]
        montecarlo_section['black_body_sampling']['end'] = \
            black_body_section[1]
        montecarlo_section['black_body_sampling']['samples'] = \
            black_body_section[2]
        ###### END of convergence section reading


        validated_config_dict['spectrum'] = parse_spectrum_list2dict(
            validated_config_dict['spectrum'])

        return cls(validated_config_dict, atom_data)


    def __init__(self, config_dict, atom_data):
        super(Configuration, self).__init__(config_dict)
        self.atom_data = atom_data
        selected_atomic_numbers = self.abundances.index
        if atom_data is not None:
            self.number_densities = (self.abundances * self.structure.mean_densities.to('g/cm^3').value)
            self.number_densities = self.number_densities.div(self.atom_data.atom_data.mass.ix[selected_atomic_numbers],
                                                              axis=0)
        else:
            logger.critical('atom_data is None, only sensible for testing the parser')







