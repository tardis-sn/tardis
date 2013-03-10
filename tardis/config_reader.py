# Module to read the rather complex config data
# Currently the configuration file is documented in
# tardis/data/example_configuration.ini

from astropy import constants, units
from ConfigParser import ConfigParser
import logging
import numpy as np
import os
import h5py
import re
import pandas as pd
from tardis import atomic
import yaml

import pdb

import pprint

logger = logging.getLogger(__name__)

data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))

#File parsers for different file formats:

density_structure_fileparser = {}


class TardisConfigError(ValueError):
    pass


def parse2quantity(quantity_string):
    value_string, unit_string = quantity_string.split()

    value = float(value_string)
    return units.Quantity(value, unit_string)


def calculate_density_after_time(densities, time_0, time_explosion):
    return densities * (time_explosion / time_0) ** -3


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


def parse_density_file_section(density_file_dict, time_explosion):
    density_file_parser = {}

    def parse_artis_density(density_file_dict, time_explosion):
        density_file = density_file_dict['name']
        for i, line in enumerate(file(density_file)):
            if i == 0:
                no_of_shells = int(line.strip())
            elif i == 1:
                time_of_model = units.Quantity(float(line.strip()), 'day').to('s')
            elif i == 2:
                break

        velocities, mean_densities_0 = np.recfromtxt(density_file, skip_header=2, usecols=(1, 2), unpack=True)
        #converting densities from log(g/cm^3) to g/cm^3 and stretching it to the current ti
        velocities = units.Quantity(np.append([0], velocities), 'km/s').to('cm/s')
        mean_densities_0 = units.Quantity(10 ** mean_densities_0, 'g/cm^3')

        mean_densities = calculate_density_after_time(mean_densities_0, time_of_model.value, time_explosion)


        #Verifying information
        if len(mean_densities) == no_of_shells:
            logger.debug('Verified ARTIS file %s (no_of_shells=length of dataset)', density_file)
        else:
            raise TardisConfigError(
                'Error in ARTIS file %s - Number of shells not the same as dataset length' % density_file)

        min_shell = 0
        max_shell = no_of_shells

        v_inner = velocities[:-1]
        v_outer = velocities[1:]

        volumes = (4 * np.pi / 3) * (time_of_model ** 3) * ( v_outer ** 3 - v_inner ** 3)
        masses = (volumes * mean_densities_0 / constants.M_sun).to(1)

        logger.info('Read ARTIS configuration file %s - found %d zones with total mass %g Msun', density_file,
                    no_of_shells, sum(masses.value))

        if 'v_lowest' in density_file_dict:
            v_lowest = parse2quantity(density_file_dict['v_lowest']).to('cm/s').value
            min_shell = v_inner.value.searchsorted(v_lowest)
        else:
            min_shell = 0

        if 'v_highest' in density_file_dict:
            v_highest = parse2quantity(density_file_dict['v_highest']).to('cm/s').value
            max_shell = v_outer.value.searchsorted(v_highest)
        else:
            max_shell = no_of_shells

        v_inner = v_inner[min_shell:max_shell]
        v_outer = v_outer[min_shell:max_shell]
        mean_densities = mean_densities[min_shell:max_shell]

        return v_inner.value, v_outer.value, mean_densities.value, min_shell, max_shell

    density_file_parser['artis'] = parse_artis_density

    try:
        parser = density_file_parser[density_file_dict['type']]
    except KeyError:
        raise TardisConfigError('In abundance file section only types %s are allowed (supplied %s) ' %
                                (density_file_parser.keys(), density_file_dict['type']))

    return parser(density_file_dict, time_explosion)


def parse_velocity_section(velocity_dict, no_of_shells):
    velocity_parser = {}

    def parse_linear_velocity(velocity_dict, no_of_shells):
        v_inner = parse2quantity(velocity_dict['v_inner']).to('cm/s').value
        v_outer = parse2quantity(velocity_dict['v_outer']).to('cm/s').value
        velocities = np.linspace(v_inner, v_outer, no_of_shells + 1)
        return velocities[:-1], velocities[1:]

    velocity_parser['linear'] = parse_linear_velocity

    try:
        parser = velocity_parser[velocity_dict['type']]
    except KeyError:
        raise TardisConfigError('In velocity section only types %s are allowed (supplied %s) ' %
                                (velocity_parser.keys(), velocity_dict['type']))
    return parser(velocity_dict, no_of_shells)


def parse_density_section(density_dict, no_of_shells, v_inner, v_outer, time_explosion):
    density_parser = {}


    #Parse density uniform
    def parse_uniform(density_dict, no_of_shells, v_inner, v_outer, time_explosion):
        return np.ones(no_of_shells) * parse2quantity(density_dict['value']).to('g cm^-3').value

    density_parser['uniform'] = parse_uniform

    #Parse density branch85 w7
    def parse_branch85(density_dict, no_of_shells, v_inner, v_outer, time_explosion):

        time_0 = density_dict.pop('time_0', 19.9999584)
        if isinstance(time_0, basestring):
            time_0 = parse2quantity(time_0).to('s').value
        else:
            logger.debug('time_0 not supplied for density branch85 - using sensible default %g', time_0)

        density_coefficient = density_dict.pop('density_coefficient', None)
        if density_coefficient is None:
            density_coefficient = 3e29
            logger.debug('density_coefficient not supplied for density type branch85 - using sensible default %g',
                         density_coefficient)

        velocities = 0.5 * (v_inner + v_outer)
        densities = density_coefficient * (velocities * 1e-5) ** -7

        densities = calculate_density_after_time(densities, time_0, time_explosion)

        return densities

    density_parser['branch85_w7'] = parse_branch85

    def parse_exponential(density_dict, no_of_shells, v_inner, v_outer, time_explosion):
        time_0 = density_dict.pop('time_0', 19.9999584)
        if isinstance(time_0, basestring):
            time_0 = parse2quantity(time_0).to('s').value
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
        raise TardisConfigError('In density section only types %s are allowed (supplied %s) ' %
                                (density_parser.keys(), density_dict['type']))
    return parser(density_dict, no_of_shells, v_inner, v_outer, time_explosion)


def parse_abundance_file_section(abundance_file_dict, abundances, min_shell, max_shell):
    abundance_file_parser = {}

    def parse_artis(abundance_file_dict, abundances, min_shell, max_shell):
        fname = abundance_file_dict['name']
        max_atom = 30
        logger.info("Parsing ARTIS Abundance section from shell %d to %d", min_shell, max_shell)
        abundances.values[:, :max_atom] = np.loadtxt(fname)[min_shell:max_shell, 1:]
        return abundances

    abundance_file_parser['artis'] = parse_artis

    try:
        parser = abundance_file_parser[abundance_file_dict['type']]
    except KeyError:
        raise TardisConfigError('In abundance file section only types %s are allowed (supplied %s) ' %
                                (abundance_file_parser.keys(), abundance_file_dict['type']))

    return parser(abundance_file_dict, abundances, min_shell, max_shell)


def reformat_element_symbol(element_symbol):
    """
    Reformat the string so the first letter is uppercase and all subsequent letters lowercase

    Parameters
    ----------
        element_symbol: str

    Returns
    -------
        reformated element symbol
    """

    #Reformating element to appropriate list
    element_str = list(element_symbol)
    element_str[0] = element_str[0].upper()
    element_str[1:] = [item.lower() for item in element_str[1:]]
    return ''.join(element_str)


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

    logger.info("Choosing uniform abundance set 'lucy99':\n %s",
                pd.DataFrame(lucy99.__array__()))

    return dict(zip(lucy99.dtype.names, lucy99[0]))


class TardisConfiguration(object):
    """
    Tardis configuration class
    """

    @classmethod
    def from_ini(cls, fname, args=None):
        print "This function is being deprecated and replaced by from_yaml classmethod"
        config_parse_object = ConfigParser()
        config_parse_object.read(fname)
        general_dict = dict(config_parse_object.items('general'))
        abundance_dict = dict(config_parse_object.items('abundances'))

        config_object = cls()
        config_object.parse_general_section(general_dict)
        config_object.parse_abundance_section(abundance_dict)
        return config_object

    @classmethod
    def from_yaml(cls, fname, args=None):
        """
        Reading in from a YAML file and commandline args. Preferring commandline args when given

        :param cls:
        :param fname:
        :param args:
        :return:
        """
        try:
            yaml_dict = yaml.load(file(fname))
        except IOError as e:
            logger.critical('No config file named: %s', fname)
            raise e
        if yaml_dict['config_type'] not in ['simple1d']:
            raise TardisConfigError('Only config_type=simple1d allowed at the moment.')

        config_dict = {}

        #First let's see if we can find an atom_db anywhere:
        if args is not None and args.atom_data is not None:
            atom_data_fname = args.atom_data
            if 'atom_data' in yaml_dict.keys():
                logger.warn('Ignoring atom_data given in config file (%s)', yaml_dict['atom_data'])
        elif 'atom_data' in yaml_dict.keys():
            atom_data_fname = yaml_dict['atom_data']
        else:
            raise TardisConfigError('No atom_data key found in config or command line')

        logger.info('Reading Atomic Data from %s', atom_data_fname)
        atom_data = atomic.AtomData.from_hdf5(atom_data_fname)
        config_dict['atom_data'] = atom_data
        #Next finding the time of explosion

        try:
            time_explosion = parse2quantity(yaml_dict['time_explosion']).to('s')
        except AttributeError as ae:
            logger.critical(str(ae))
            raise ae

        config_dict['time_explosion'] = time_explosion

        if 'log_lsun' in yaml_dict['luminosity']:
            luminosity_value, luminosity_unit = yaml_dict['luminosity'].split()
            config_dict['luminosity'] = 10 ** (float(luminosity_value) + np.log10(constants.L_sun.cgs.value))
        else:
            config_dict['luminosity'] = parse2quantity(yaml_dict['luminosity'])


        #Trying to figure out the structure (number of shells)
        structure_dict = yaml_dict['model'].pop('structure')

        #first let's try to see if there's a file keyword
        if 'file' in structure_dict.keys():
            density_file_section = structure_dict.pop('file')
            v_inner, v_outer, mean_densities, min_shell, max_shell = parse_density_file_section(density_file_section,
                                                                                                time_explosion)

            no_of_shells = len(v_inner)
            if structure_dict != {}:
                logger.warn(
                    'Accepted file for structure (density/velocity) structure ignoring all other arguments: \n%s\n',
                    pprint.pformat(structure_dict, indent=4))
        else:
            #requiring all keys: no_of_shells, velocity, density
            if not all([item in structure_dict.keys() for item in ('no_of_shells', 'velocity', 'density')]):
                raise TardisConfigError(
                    'If file-section is not given to structure-section, one needs to provide all: no_of_shells, velocity, density')

            no_of_shells = structure_dict['no_of_shells']

            v_inner, v_outer = parse_velocity_section(structure_dict['velocity'], no_of_shells)
            mean_densities = parse_density_section(structure_dict['density'], no_of_shells, v_inner, v_outer,
                                                   time_explosion)

        config_dict['v_inner'] = v_inner
        config_dict['v_outer'] = v_outer
        config_dict['mean_densities'] = mean_densities
        config_dict['no_of_shells'] = no_of_shells

        #Now that the structure section is parsed we move on to the abundances

        abundances_dict = yaml_dict['model']['abundances'].copy()
        #TODO: columns are now until Z=120
        species_pattern = re.compile('\s*([a-zA-Z]*)(\d*)\s*')
        abundances = pd.DataFrame(columns=np.arange(1, 120), index=pd.Index(np.arange(no_of_shells), name='shells'))

        if 'file' in abundances_dict.keys():
            abundance_file_dict = abundances_dict.pop('file')
            parse_abundance_file_section(abundance_file_dict, abundances, min_shell, max_shell)

        if 'abundance_set' in abundances_dict.keys():
            abundance_set_dict = abundances_dict.pop('abundance_set')
            print "abundance set not implemented currently"
        #            abundance_set = abundance_dict.pop('abundance_set', None)
        #            if abundance_set == 'lucy99':
        #                abundances = read_lucy99_abundances()
        #            elif abundance_set is not None:
        #                raise ValueError("Currently only 'lucy99' abundance_set supported")

        nlte_species = []
        if 'nlte_species' in abundances_dict.keys():
            nlte_species_list = abundances_dict.pop('nlte_species')
            for species_symbol in nlte_species_list:
                species_match = species_pattern.match(species_symbol)
                if species_match is None:
                    raise ValueError(
                        "'nlte_species' element %s could not be matched to a valid species notation (e.g. Ca2)")
                species_element, species_ion = species_match.groups()
                species_element = reformat_element_symbol(species_element)
                if species_element not in atom_data.symbol2atomic_number:
                    raise ValueError("Element provided in NLTE species %s unknown" % species_element)
                nlte_species.append((atom_data.symbol2atomic_number[species_element], int(species_ion) - 1))

        for element in abundances_dict:
            element_symbol = reformat_element_symbol(element)
            if element_symbol not in atom_data.symbol2atomic_number:
                raise ValueError('Element %s provided in config unknown' % element_symbol)

            z = atom_data.symbol2atomic_number[element_symbol]

            abundances[z] = float(abundances_dict[element])

        config_dict['abundances'] = abundances
        config_dict['nlte_species'] = nlte_species


        ########### DOING PLASMA SECTION ###############

        plasma_section = yaml_dict.pop('plasma')

        config_dict['initial_t_rad'] = parse2quantity(plasma_section['initial_t_rad']).to('K').value
        config_dict['initial_t_inner'] = parse2quantity(plasma_section['initial_t_inner']).to('K').value

        if plasma_section['plasma_type'] not in ('nebular', 'lte'):
            raise TardisConfigError('plasma_type only allowed to be "nebular" or "lte"')
        config_dict['plasma_type'] = plasma_section['plasma_type']

        if plasma_section['radiative_rates_type'] not in ('nebular', 'lte', 'detailed'):
            raise TardisConfigError('radiative_rates_types must be either "nebular", "lte", or "detailed"')
        config_dict['radiative_rates_type'] = plasma_section['radiative_rates_type']

        if plasma_section['line_interaction_type'] not in ('scatter', 'downbranch', 'macroatom'):
            raise TardisConfigError('radiative_rates_types must be either "scatter", "downbranch", or "macroatom"')
        config_dict['line_interaction_type'] = plasma_section['line_interaction_type']

        config_dict.update(yaml_dict.pop('montecarlo', {}))

        disable_electron_scattering = plasma_section['disable_electron_scattering']

        if disable_electron_scattering is False:
            logger.info("Electron scattering switched on")
            config_dict['sigma_thomson'] = None
        else:
            logger.warn('Disabling electron scattering - this is not physical')
            config_dict['sigma_thomson'] = 1e-200


            ##### spectrum section ######
        spectrum_section = yaml_dict.pop('spectrum')
        spectrum_start = parse2quantity(spectrum_section['start']).to('angstrom', units.spectral())
        spectrum_end = parse2quantity(spectrum_section['end']).to('angstrom', units.spectral())
        spectrum_bins = int(spectrum_section['bins'])

        if spectrum_end > spectrum_start:
            logger.debug('Converted spectrum start/end to angstrom %.4g %.4g', spectrum_start, spectrum_end)
            spectrum_start = spectrum_start
            spectrum_end = spectrum_end

        else:
            logger.warn('Spectrum Start > Spectrum End in wavelength space - flipped them')

            logger.debug('Converted spectrum start/end to angstrom %.4g %.4g', spectrum_end, spectrum_start)
            tmp = spectrum_start
            spectrum_start = spectrum_end
            spectrum_end = tmp
        config_dict['spectrum_start'] = spectrum_start
        config_dict['spectrum_end'] = spectrum_end
        config_dict['spectrum_bins'] = spectrum_bins

        config_dict['spectrum_start_nu'] = spectrum_end.to('Hz', units.spectral())
        config_dict['spectrum_end_nu'] = spectrum_start.to('Hz', units.spectral())

        sn_distance = spectrum_section.pop('sn_distance', None)

        if sn_distance is not None:
            if sn_distance.strip().lower() == 'lum_density':
                logger.info('Luminosity density requested  - setting distance to sqrt(1/(4*pi))')
                config_dict['lum_density'] = True
                config_dict['sn_distance'] = np.sqrt(1 / (4 * np.pi))
            else:
                config_dict['sn_distance'] = parse2quantity(sn_distance).to('cm').value
                config_dict['lum_density'] = False
        else:
            config_dict['sn_distance'] = None

        return cls(config_dict)


    def __init__(self, config_dict):

        for key in config_dict:
            setattr(self, key, config_dict[key])

        self.number_densities = self.calculate_number_densities()


    def calculate_number_densities(self):
        abundances = self.abundances
        for atomic_number in abundances:

            if all(abundances[atomic_number].isnull()):
                del abundances[atomic_number]
                continue
            else:
                abundances[abundances[atomic_number].isnull()] == 0.0

        #normalizing
        abundances = abundances.divide(abundances.sum(axis=1), axis=0)
        atom_mass = self.atom_data.atom_data.ix[abundances.columns].mass
        number_densities = (abundances.mul(self.mean_densities, axis=0)).divide(atom_mass)

        return number_densities





