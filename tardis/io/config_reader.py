# Module to read the rather complex config data
# Currently the configuration file is documented in
# tardis/data/example_configuration.ini

import logging
import os
import pprint
import copy

from astropy import constants, units as u
import astropy.utils
import numpy as np
import pandas as pd
import yaml
from model_reader import read_density_file, calculate_density_after_time, read_abundances_file

from tardis import atomic
from tardis.util import species_string_to_tuple, parse_quantity, element_symbol2atomic_number

pp = pprint.PrettyPrinter(indent=4)

logger = logging.getLogger(__name__)

data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../data'))

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

    Parameters:
    -----------

    quantity_linspace_dictionary: ~dict

    add_one: boolean, default: True

    Returns:
    --------

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



def calculate_exponential_densities(velocities, velocity_0, rho_0, exponent):
    """

    This function computes a descret exponential density profile.
    :math:`\\rho = \\rho_0 \\times \\left( \\frac{v_0}{v} \\right)^n`

    Parameters
    ----------

    velocities : Array like list
                velocities in km/s

    velocity_0 : ~float
        Velocity at the inner boundary


    rho_0 : ~float
        density at velocity_0

    exponent : ~float
        exponent used in the powerlaw

    Returns
    -------

    Array like density structure

    """
    densities = rho_0 * np.power((velocity_0 / velocities), exponent)
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
        abundances = pd.DataFrame(np.loadtxt(abundances_fname)[min_shell:max_shell, 1:].transpose(), index=np.arange(1, 31))

        ni_stable = abundances.ix[28] - artis_model['ni56_fraction']
        co_stable = abundances.ix[27] - artis_model['co56_fraction']
        fe_stable = abundances.ix[26] - artis_model['fe52_fraction']
        mn_stable = abundances.ix[25] - 0.0
        cr_stable = abundances.ix[24] - artis_model['cr48_fraction']
        v_stable = abundances.ix[23] - 0.0
        ti_stable = abundances.ix[22] - 0.0


        abundances.ix[28] = ni_stable
        abundances.ix[28] += artis_model['ni56_fraction'] * np.exp(-(time_explosion* inv_ni56_efolding_time).to(1).value)

        abundances.ix[27] = co_stable
        abundances.ix[27] += artis_model['co56_fraction'] * np.exp(-(time_explosion* inv_co56_efolding_time).to(1).value)
        abundances.ix[27] += (inv_ni56_efolding_time * artis_model['ni56_fraction'] /
                              (inv_ni56_efolding_time - inv_co56_efolding_time)) * \
                             (np.exp(-(inv_co56_efolding_time * time_explosion).to(1).value) - np.exp(-(inv_ni56_efolding_time * time_explosion).to(1).value))

        abundances.ix[26] = fe_stable
        abundances.ix[26] += artis_model['fe52_fraction'] * np.exp(-(time_explosion * inv_fe52_efolding_time).to(1).value)
        abundances.ix[26] += ((artis_model['co56_fraction'] * inv_ni56_efolding_time
                               - artis_model['co56_fraction'] * inv_co56_efolding_time
                               + artis_model['ni56_fraction'] * inv_ni56_efolding_time
                               - artis_model['ni56_fraction'] * inv_co56_efolding_time
                               - artis_model['co56_fraction'] * inv_ni56_efolding_time * np.exp(-(inv_co56_efolding_time * time_explosion).to(1).value)
                               + artis_model['co56_fraction'] * inv_co56_efolding_time * np.exp(-(inv_co56_efolding_time * time_explosion).to(1).value)
                               - artis_model['ni56_fraction'] * inv_ni56_efolding_time * np.exp(-(inv_co56_efolding_time * time_explosion).to(1).value)
                               + artis_model['ni56_fraction'] * inv_co56_efolding_time * np.exp(-(inv_ni56_efolding_time * time_explosion).to(1).value))
        / (inv_ni56_efolding_time - inv_co56_efolding_time))


        abundances.ix[25] = mn_stable
        abundances.ix[25] += (inv_fe52_efolding_time * artis_model['fe52_fraction'] /
                              (inv_fe52_efolding_time - inv_mn52_efolding_time)) * \
                             (np.exp(-(inv_mn52_efolding_time * time_explosion).to(1).value) - np.exp(-(inv_fe52_efolding_time * time_explosion).to(1).value))

        abundances.ix[24] = cr_stable
        abundances.ix[24] += artis_model['cr48_fraction'] * np.exp(-(time_explosion* inv_cr48_efolding_time).to(1).value)
        abundances.ix[24] += ((artis_model['fe52_fraction'] * inv_fe52_efolding_time
                               - artis_model['fe52_fraction'] * inv_mn52_efolding_time
                               - artis_model['fe52_fraction'] * inv_fe52_efolding_time * np.exp(-(inv_mn52_efolding_time * time_explosion).to(1).value)
                               + artis_model['fe52_fraction'] * inv_mn52_efolding_time * np.exp(-(inv_fe52_efolding_time * time_explosion).to(1).value))
        / (inv_fe52_efolding_time - inv_mn52_efolding_time))

        abundances.ix[23] = v_stable
        abundances.ix[23] += (inv_cr48_efolding_time * artis_model['cr48_fraction'] /
                              (inv_cr48_efolding_time - inv_v48_efolding_time)) * \
                             (np.exp(-(inv_v48_efolding_time * time_explosion).to(1).value) - np.exp(-(inv_cr48_efolding_time * time_explosion).to(1).value))

        abundances.ix[22] = ti_stable
        abundances.ix[22] += ((artis_model['cr48_fraction'] * inv_cr48_efolding_time
                               - artis_model['cr48_fraction'] * inv_v48_efolding_time
                               - artis_model['cr48_fraction'] * inv_cr48_efolding_time * np.exp(-(inv_v48_efolding_time * time_explosion).to(1).value)
                               + artis_model['cr48_fraction'] * inv_v48_efolding_time * np.exp(-(inv_cr48_efolding_time * time_explosion).to(1).value))
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
            mean_densities = np.empty(no_of_shells*split_shells) * old_mean_densities.unit
            new_abundance_data = np.empty((abundances.values.shape[0], no_of_shells * split_shells))
            for i in xrange(split_shells):
                mean_densities[i::split_shells] = old_mean_densities
                new_abundance_data[:,i::split_shells] = abundances.values

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
        return parse_quantity(density_dict['value']).to('g cm^-3') * np.ones(no_of_shells)

    density_parser['uniform'] = parse_uniform

    #Parse density branch85 w7
    def parse_branch85(density_dict, v_inner, v_outer, time_explosion):

        time_0 = density_dict.pop('time_0', 19.9999584)
        if isinstance(time_0, basestring):
            time_0 = parse_quantity(time_0).to('s')
        else:
            time_0 *= u.s
            logger.debug('time_0 not supplied for density branch85 - using sensible default %g', time_0)

        density_coefficient = density_dict.pop('density_coefficient', None)
        if density_coefficient is None:
            density_coefficient = 3e29 * u.Unit('g/cm^3')
            logger.debug('density_coefficient not supplied for density type branch85 - using sensible default %g',
                         density_coefficient)
        else:
            density_coefficient = parse_quantity(density_coefficient)

        velocities = 0.5 * (v_inner + v_outer)
        densities = density_coefficient * (velocities.value * 1e-5) ** -7

        densities = calculate_density_after_time(densities, time_0, time_explosion)

        return densities

    density_parser['branch85_w7'] = parse_branch85

    def parse_exponential(density_dict, v_inner, v_outer, time_explosion):
        time_0 = density_dict.pop('time_0', 19.9999584)
        if isinstance(time_0, basestring):
            time_0 = parse_quantity(time_0).to('s').value
        else:
            logger.debug('time_0 not supplied for density branch85 - using sensible default %g', time_0)
        try:
            rho_0 = float(density_dict.pop('rho_0'))
        except KeyError:
            rho_0 = 1e-2
            logger.warning('rho_o was not given in the config! Using %g', rho_0)
        try:
            exponent = density_dict.pop('exponent')
        except KeyError:
            exponent = 2
            logger.warning('exponent was not given in the config file! Using %f', exponent)

        velocities = 0.5 * (v_inner + v_outer)
        densities = calculate_exponential_densities(velocities, v_inner[0], rho_0, exponent)

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
        config_dict['luminosity_requested'] = 10 ** (float(luminosity_value) + np.log10(constants.L_sun.cgs.value)) * u.erg / u.s
    else:
        config_dict['luminosity_requested'] = (float(luminosity_value) * u.Unit(luminosity_unit)).to('erg/s')

    config_dict['time_explosion'] = parse_quantity(supernova_dict['time_explosion']).to('s')

    if 'distance' in supernova_dict:
        config_dict['distance'] = parse_quantity(supernova_dict['distance'])
    else:
        config_dict['distance'] = None

    if 'luminosity_wavelength_start' in supernova_dict:
        config_dict['luminosity_nu_end'] = parse_quantity(supernova_dict['luminosity_wavelength_start']).\
            to('Hz', u.spectral())
    else:
        config_dict['luminosity_nu_end'] = np.inf * u.Hz

    if 'luminosity_wavelength_end' in supernova_dict:
        config_dict['luminosity_nu_start'] = parse_quantity(supernova_dict['luminosity_wavelength_end']).\
            to('Hz', u.spectral())
    else:
        config_dict['luminosity_nu_start'] = 0.0 * u.Hz

    return config_dict



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





class TARDISConfigurationNameSpace(object):

    def __init__(self, config_dict):
        self.config_dict = config_dict


    def __getattr__(self, item):
        if item in self.config_dict:
            config_item = self.config_dict[item]
            if isinstance(config_item, dict):
                setattr(self, item, TARDISConfigurationNameSpace(config_item))
                return getattr(self, item)
            else:
                return self.config_dict[item]
        else:
            return super(TARDISConfigurationNameSpace, self).__getattribute__(item)

    def __getitem__(self, item):
        return self.config_dict.__getitem__(item)

    def get(self, k, d=None):
        return self.config_dict.get(k, d)


    def __repr__(self):
        return pp.pformat(self.config_dict)

    def __dir__(self):
        return self.__dict__.keys() + self.config_dict.keys()


class TARDISConfiguration(TARDISConfigurationNameSpace):
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
    def from_config_dict(cls, raw_dict, atom_data=None, test_parser=False):
        """
        Reading in from a YAML file and commandline args. Preferring commandline args when given

        Parameters
        ----------

        fname : filename for the yaml file

        args : namespace object
            Not implemented Yet

        Returns
        -------

        `tardis.config_reader.TARDISConfiguration`

        """

        config_dict = {}
        raw_dict = copy.deepcopy(raw_dict)

        #First let's see if we can find an atom_db anywhere:
        if test_parser:
          atom_data = None
        elif 'atom_data' in raw_dict.keys():
            atom_data_fname = raw_dict['atom_data']
            config_dict['atom_data_fname'] = atom_data_fname
        else:
            raise ConfigurationError('No atom_data key found in config or command line')



        if atom_data is None and not test_parser:
            logger.info('Reading Atomic Data from %s', atom_data_fname)
            atom_data = atomic.AtomData.from_hdf5(atom_data_fname)
        else:
            atom_data = atom_data



        #Parsing supernova dictionary
        config_dict['supernova'] = parse_supernova_section(raw_dict['supernova'])

        #Parsing the model section
        model_section = raw_dict.pop('model')
        v_inner = None
        v_outer = None
        mean_densities = None
        abundances = None


        if 'file' in model_section:
            v_inner, v_outer, mean_densities, abundances = parse_model_file_section(model_section.pop('file'),
                                                                                    config_dict['supernova']['time_explosion'])
            no_of_shells = len(v_inner)

        structure_config_dict = {}

        if 'structure' in model_section:
        #Trying to figure out the structure (number of shells)

            structure_section = model_section.pop('structure')
            inner_boundary_index, outer_boundary_index = None, None
            try:
                structure_section_type = structure_section['type']
            except KeyError:
                raise ConfigurationError('Structure section requires "type" keyword')


            if structure_section_type == 'specific':
                velocities = parse_quantity_linspace(structure_section['velocity']).to('cm/s')
                v_inner, v_outer = velocities[:-1], velocities[1:]

                mean_densities = parse_density_section(structure_section['density'], v_inner, v_outer,
                                                       config_dict['supernova']['time_explosion'])

            elif structure_section_type == 'file':
                v_inner_boundary, v_outer_boundary = structure_section.get('v_inner_boundary', 0 * u.km/u.s), \
                                                     structure_section.get('v_outer_boundary', np.inf * u.km/u.s)

                if not hasattr(v_inner_boundary, 'unit'):
                    v_inner_boundary = parse_quantity(v_inner_boundary)

                if not hasattr(v_outer_boundary, 'unit'):
                    v_outer_boundary = parse_quantity(v_outer_boundary)

                v_inner, v_outer, mean_densities, inner_boundary_index, outer_boundary_index =\
                    read_density_file(structure_section['filename'], structure_section['filetype'],
                                      config_dict['supernova']['time_explosion'], v_inner_boundary, v_outer_boundary)
        else:
            raise ConfigurationError('structure section required in configuration file')


        r_inner = config_dict['supernova']['time_explosion'] * v_inner
        r_outer = config_dict['supernova']['time_explosion'] * v_outer
        r_middle = 0.5 * (r_inner + r_outer)

        structure_config_dict['v_inner'] = v_inner
        structure_config_dict['v_outer'] = v_outer
        structure_config_dict['mean_densities'] = mean_densities
        no_of_shells = len(v_inner)
        structure_config_dict['no_of_shells'] = no_of_shells
        structure_config_dict['r_inner'] = r_inner
        structure_config_dict['r_outer'] = r_outer
        structure_config_dict['r_middle'] = r_middle
        structure_config_dict['volumes'] = (4. / 3) * np.pi * (r_outer ** 3 - r_inner ** 3)




        config_dict['structure'] = structure_config_dict
        #Now that the structure section is parsed we move on to the abundances



        abundances_section  = model_section.pop('abundances')
        abundances_type = abundances_section.pop('type')

        if abundances_type == 'uniform':
            abundances = pd.DataFrame(columns=np.arange(no_of_shells),
                  index=pd.Index(np.arange(1, 120), name='atomic_number'), dtype=np.float64)

            for element_symbol_string in abundances_section:

                z = element_symbol2atomic_number(element_symbol_string)
                abundances.ix[z] = float(abundances_section[element_symbol_string])

        elif abundances_type == 'file':
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

        config_dict['abundances'] = abundances



        ########### DOING PLASMA SECTION ###############

        plasma_section = raw_dict.pop('plasma')
        plasma_config_dict = {}

        if plasma_section['ionization'] not in ('nebular', 'lte'):
            raise ConfigurationError('plasma_type only allowed to be "nebular" or "lte"')
        plasma_config_dict['ionization'] = plasma_section['ionization']


        if plasma_section['excitation'] not in ('dilute-lte', 'lte'):
            raise ConfigurationError('plasma_type only allowed to be "nebular" or "lte"')
        plasma_config_dict['excitation'] = plasma_section['excitation']

        if plasma_section['radiative_rates_type'] not in ('dilute-blackbody', 'detailed'):
            raise ConfigurationError('radiative_rates_types must be either "dilute-blackbody" or "detailed"')
        plasma_config_dict['radiative_rates_type'] = plasma_section['radiative_rates_type']

        if plasma_section['line_interaction_type'] not in ('scatter', 'downbranch', 'macroatom'):
            raise ConfigurationError('radiative_rates_types must be either "scatter", "downbranch", or "macroatom"')
        plasma_config_dict['line_interaction_type'] = plasma_section['line_interaction_type']

        if 'w_epsilon' in plasma_section:
            plasma_config_dict['w_epsilon'] = plasma_section['w_epsilon']
        else:
            logger.warn('"w_epsilon" not specified in plasma section - setting it to 1e-10')
            plasma_config_dict['w_epsilon'] = 1e-10

        if 'delta_treatment' in plasma_section:
            plasma_config_dict['delta_treatment'] = plasma_section['delta_treatment']
        else:
            logger.warn('"delta_treatment" not specified in plasma section - defaulting to None')
            plasma_config_dict['delta_treatment'] = None

        if 'initial_t_inner' in plasma_section:
            plasma_config_dict['t_inner'] = parse_quantity(plasma_section['initial_t_inner']).to('K')
        else:
            plasma_config_dict['t_inner'] = (((config_dict['supernova']['luminosity_requested'] / \
                                            (4 * np.pi * r_inner[0]**2 * constants.sigma_sb))**.5)**.5).to('K')
            logger.info('"initial_t_inner" is not specified in the plasma section - '
                        'initializing to %s with given luminosity', plasma_config_dict['t_inner'])

        if 'initial_t_rads' in plasma_section:
            if isinstance('initial_t_rads', basestring):
                    uniform_t_rads = parse_quantity(plasma_section['initial_t_rads'])
                    plasma_config_dict['t_rads'] = u.Quantity(np.ones(no_of_shells) * uniform_t_rads.value, u.K)

            elif astropy.utils.isiterable(plasma_section['initial_t_rads']):
                assert len(plasma_section['initial_t_rads']) == no_of_shells
                plasma_config_dict['t_rads'] = u.Quantity(plasma_section['initial_t_rads'], u.K)
        else:
            logger.info('No "initial_t_rads" specified - initializing with 10000 K')

            plasma_config_dict['t_rads'] =  u.Quantity(np.ones(no_of_shells) * 10000., u.K)

        ##### NLTE subsection of Plasma start
        nlte_config_dict = {}
        nlte_species = []
        if 'nlte' in plasma_section:
            nlte_section = plasma_section['nlte']
            if 'species' in nlte_section:
                nlte_species_list = nlte_section.pop('species')
                for species_string in nlte_species_list:
                    nlte_species.append(species_string_to_tuple(species_string))

                nlte_config_dict['species'] = nlte_species
                nlte_config_dict['species_string'] = nlte_species_list
                nlte_config_dict.update(nlte_section)

                if 'coronal_approximation' not in nlte_section:
                    logger.debug('NLTE "coronal_approximation" not specified in NLTE section - defaulting to False')
                    nlte_config_dict['coronal_approximation'] = False

                if 'classical_nebular' not in nlte_section:
                    logger.debug('NLTE "classical_nebular" not specified in NLTE section - defaulting to False')
                    nlte_config_dict['classical_nebular'] = False


            elif nlte_section: #checks that the dictionary is not empty
                logger.warn('No "species" given - ignoring other NLTE options given:\n%s',
                            pp.pformat(nlte_section))

        if not nlte_config_dict:
            nlte_config_dict['species'] = []

        plasma_config_dict['nlte'] = nlte_config_dict



        #^^^^^^^ NLTE subsection of Plasma end

        config_dict['plasma'] = plasma_config_dict


        #^^^^^^^^^^^^^^ End of Plasma Section

        ##### Monte Carlo Section

        montecarlo_section = raw_dict.pop('montecarlo')
        montecarlo_config_dict = {}

        #PARSING convergence section
        convergence_variables = ['t_inner', 't_rad', 'w']
        convergence_config_dict = {}
        if 'convergence_strategy' in montecarlo_section:

            convergence_section = montecarlo_section.pop('convergence_strategy')
            if 'lock_t_inner_cycles' in convergence_section:
                lock_t_inner_cycles = convergence_section['lock_t_inner_cycles']
                logger.info('lock_t_inner_cycles set to %d cycles', lock_t_inner_cycles)
            else:
                lock_t_inner_cycles = None

            if 't_inner_update_exponent' in convergence_section:
                t_inner_update_exponent = convergence_section['t_inner_update_exponent']
                logger.info('t_inner update exponent set to %g', t_inner_update_exponent)
            else:
                t_inner_update_exponent = None

            if convergence_section['type'] == 'damped':
                convergence_config_dict['type'] == 'damped'
                global_damping_constant = convergence_section['damping_constant']

                for convergence_variable in convergence_variables:
                    convergence_parameter_name = convergence_variable
                    current_convergence_parameters = {}
                    convergence_config_dict[convergence_parameter_name] = current_convergence_parameters

                    if convergence_variable in convergence_section:
                        current_convergence_parameters['damping_constant'] \
                            = convergence_section[convergence_variable]['damping_constant']
                    else:
                        current_convergence_parameters['damping_constant'] = global_damping_constant

            elif convergence_section['type'] == 'specific':

                convergence_config_dict['type'] = 'specific'

                global_convergence_parameters = {}
                global_convergence_parameters['damping_constant'] = convergence_section['damping_constant']
                global_convergence_parameters['threshold'] = convergence_section['threshold']

                global_convergence_parameters['fraction'] = convergence_section['fraction']

                for convergence_variable in convergence_variables:
                    convergence_parameter_name = convergence_variable
                    current_convergence_parameters = {}

                    convergence_config_dict[convergence_parameter_name] = current_convergence_parameters
                    if convergence_variable in convergence_section:
                        for param in global_convergence_parameters.keys():
                            if param == 'fraction' and convergence_variable == 't_inner':
                                continue
                            if param in convergence_section[convergence_variable]:
                                current_convergence_parameters[param] = convergence_section[convergence_variable][param]
                            else:
                                current_convergence_parameters[param] = global_convergence_parameters[param]
                    else:
                        convergence_config_dict[convergence_parameter_name] = global_convergence_parameters.copy()

                global_convergence_parameters['hold'] = convergence_section['hold']
                convergence_config_dict['global_convergence_parameters'] = global_convergence_parameters

            else:
                raise ValueError("convergence criteria unclear %s", convergence_section['type'])



        else:
            lock_t_inner_cycles = None
            t_inner_update_exponent = None
            logger.warning('No convergence criteria selected - just damping by 0.5 for w, t_rad and t_inner')
            convergence_config_dict['type'] = 'damped'
            for convergence_variable in convergence_variables:
                convergence_parameter_name = convergence_variable
                convergence_config_dict[convergence_parameter_name] = dict(damping_constant=0.5)
        if lock_t_inner_cycles is None:
            logger.warning('t_inner update lock cycles not set - defaulting to 1')
            lock_t_inner_cycles = 1
        if t_inner_update_exponent is None:
            logger.warning('t_inner update exponent not set - defaulting to -0.5')
            t_inner_update_exponent = -0.5

        convergence_config_dict['lock_t_inner_cycles'] = lock_t_inner_cycles
        convergence_config_dict['t_inner_update_exponent'] = t_inner_update_exponent


        montecarlo_config_dict['convergence'] = convergence_config_dict
        ###### END of convergence section reading

        if 'last_no_of_packets' not in montecarlo_section:
            montecarlo_section['last_no_of_packets'] = None

        if 'no_of_virtual_packets' not in montecarlo_section:
            montecarlo_section['no_of_virtual_packets'] = 0

        montecarlo_config_dict.update(montecarlo_section)

        disable_electron_scattering = plasma_section.get('disable_electron_scattering', False)

        if disable_electron_scattering is False:
            logger.info("Electron scattering switched on")
            montecarlo_config_dict['sigma_thomson'] =6.652486e-25 / (u.cm**2)
        else:
            logger.warn('Disabling electron scattering - this is not physical')
            montecarlo_config_dict['sigma_thomson'] = 1e-200 / (u.cm**2)

        montecarlo_config_dict['enable_reflective_inner_boundary'] = False
        montecarlo_config_dict['inner_boundary_albedo'] = 0.0

        if 'inner_boundary_albedo' in montecarlo_section:
            montecarlo_config_dict['inner_boundary_albedo'] = montecarlo_section['inner_boundary_albedo']
            if 'enable_reflective_inner_boundary' not in montecarlo_section:
                logger.warn('inner_boundary_albedo set, however enable_reflective_inner_boundary option not specified '
                            '- defaulting to reflective inner boundary')
                montecarlo_config_dict['enable_reflective_inner_boundary'] = True

            if 'enable_reflective_inner_boundary' in montecarlo_section:
                montecarlo_config_dict['enable_reflective_inner_boundary'] = montecarlo_section['enable_reflective_inner_boundary']
                if montecarlo_section['enable_reflective_inner_boundary'] == True and 'inner_boundary_albedo' not in montecarlo_section:
                    logger.warn('enabled reflective inner boundary, but "inner_boundary_albedo" not set - defaulting to 0.5')
                    montecarlo_config_dict['inner_boundary_albedo'] = 0.5




        if 'black_body_sampling' in montecarlo_section:
            black_body_sampling_section = montecarlo_section.pop('black_body_sampling')
            sampling_start, sampling_end = parse_spectral_bin(black_body_sampling_section['start'],
                                                                                black_body_sampling_section['stop'])
            montecarlo_config_dict['black_body_sampling']['start'] = sampling_start
            montecarlo_config_dict['black_body_sampling']['end'] = sampling_end
            montecarlo_config_dict['black_body_sampling']['samples'] = np.int64(black_body_sampling_section['num'])
        else:
            logger.warn('No "black_body_sampling" section in config file - using defaults of '
                        '50 - 200000 Angstrom (1e6 samples)')
            montecarlo_config_dict['black_body_sampling'] = {}
            montecarlo_config_dict['black_body_sampling']['start'] = 50 * u.angstrom
            montecarlo_config_dict['black_body_sampling']['end'] = 200000 * u.angstrom
            montecarlo_config_dict['black_body_sampling']['samples'] = np.int64(1e6)

        config_dict['montecarlo'] = montecarlo_config_dict
        ##### End of MonteCarlo section






        ##### spectrum section ######
        spectrum_section = raw_dict.pop('spectrum')
        spectrum_config_dict = {}
        spectrum_frequency = parse_quantity_linspace(spectrum_section).to('Hz', u.spectral())

        if spectrum_frequency[0] > spectrum_frequency[1]:
            spectrum_frequency = spectrum_frequency[::-1]

        spectrum_config_dict['start'] = parse_quantity(spectrum_section['start'])
        spectrum_config_dict['end'] = parse_quantity(spectrum_section['stop'])
        spectrum_config_dict['bins'] = spectrum_section['num']


        spectrum_frequency = np.linspace(spectrum_config_dict['end'].to('Hz', u.spectral()).value,
                                                         spectrum_config_dict['start'].to('Hz', u.spectral()).value,
                                                         num=spectrum_config_dict['bins'] + 1) * u.Hz

        spectrum_config_dict['frequency'] = spectrum_frequency.to('Hz')
        config_dict['spectrum'] = spectrum_config_dict




        return cls(config_dict, atom_data)


    def __init__(self, config_dict, atom_data):
        super(TARDISConfiguration, self).__init__(config_dict)
        self.atom_data = atom_data
        selected_atomic_numbers = self.abundances.index
        if atom_data is not None:
            self.number_densities = (self.abundances * self.structure.mean_densities.to('g/cm^3').value)
            self.number_densities = self.number_densities.div(self.atom_data.atom_data.mass.ix[selected_atomic_numbers],
                                                              axis=0)
        else:
            logger.critical('atom_data is None, only sensible for testing the parser')







