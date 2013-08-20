# Module to read the rather complex config data
# Currently the configuration file is documented in
# tardis/data/example_configuration.ini

import logging
import os
import pprint

from astropy import constants, units as u
import astropy.utils
import numpy as np
import h5py
import pandas as pd
import yaml

import tardis.util
from tardis import atomic

pp = pprint.PrettyPrinter(indent=4)

logger = logging.getLogger(__name__)

data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))

#File parsers for different file formats:

density_structure_fileparser = {}


class TARDISConfigurationError(ValueError):
    pass

class TARDISMalformedQuantityError(TARDISConfigurationError):

    def __init__(self, malformed_quantity_string):
        self.malformed_quantity_string = malformed_quantity_string

    def __str__(self):
        return 'Expecting a quantity string(e.g. "5 km/s") for keyword - supplied %s' % self.malformed_quantity_string


class TARDISMalformedElementSymbolError(TARDISConfigurationError):

    def __init__(self, malformed_element_symbol):
        self.malformed_element_symbol = malformed_element_symbol

    def __str__(self):
        return 'Expecting an atomic symbol (e.g. Fe) - supplied %s' % self.malformed_element_symbol


class TARDISMalformedSpeciesError(TARDISConfigurationError):

    def __init__(self, malformed_element_symbol):
        self.malformed_element_symbol = malformed_element_symbol

    def __str__(self):
        return 'Expecting a species notation (e.g. "Si 2", "Si II", "Fe IV") - supplied %s' % self.malformed_element_symbol



def parse_quantity(quantity_string):

    if not isinstance(quantity_string, basestring):
        raise TARDISMalformedQuantityError(quantity_string)

    try:
        value_string, unit_string = quantity_string.split()
    except ValueError:
        raise TARDISMalformedQuantityError(quantity_string)

    try:
        value = float(value_string)
    except ValueError:
        raise TARDISMalformedQuantityError(quantity_string)

    try:
        q = u.Quantity(value, unit_string)
    except ValueError:
        raise TARDISMalformedQuantityError(quantity_string)

    return q

def element_symbol2atomic_number(element_string, atom_data):
    reformatted_element_string = reformat_element_symbol(element_string)
    if reformatted_element_string not in atom_data.symbol2atomic_number:
        raise TARDISMalformedElementSymbolError(element_string)
    return atom_data.symbol2atomic_number[reformatted_element_string]


def species_tuple_to_string(species_tuple, atom_data, roman_numerals=True):
    atomic_number, ion_number = species_tuple
    element_symbol = atom_data.atomic_number2symbol[atomic_number]
    if roman_numerals:
        roman_ion_number = tardis.util.int_to_roman(ion_number+1)
        return '%s %s' % (element_symbol, roman_ion_number)
    else:
        return '%s %d' % (element_symbol, ion_number)

def species_string_to_tuple(species_string, atom_data):
    try:
        element_string, ion_number_string = species_string.split()
    except ValueError:
        raise TARDISMalformedElementSymbolError(species_string)

    atomic_number = element_symbol2atomic_number(element_string, atom_data)

    try:
        ion_number = tardis.util.roman_to_int(ion_number_string.strip())
    except ValueError:
        try:
            ion_number = int(ion_number_string)
        except ValueError:
            raise TARDISMalformedSpeciesError
    if ion_number > atomic_number:
        raise TARDISConfigurationError('Species given does not exist: ion number > atomic number')

    return atomic_number, ion_number-1





def reformat_element_symbol(element_string):
    """
    Reformat the string so the first letter is uppercase and all subsequent letters lowercase

    Parameters
    ----------
        element_symbol: str

    Returns
    -------
        reformated element symbol
    """

    return element_string[0].upper() + element_string[1:].lower()

def parse_abundance_dict_to_dataframe(abundance_dict, atom_data):
    atomic_number_dict = dict([(element_symbol2atomic_number(symbol, atom_data), abundance_dict[symbol])
                                   for symbol in abundance_dict])
    atomic_numbers = sorted(atomic_number_dict.keys())

    abundances = pd.Series([atomic_number_dict[z] for z in atomic_numbers], index=atomic_numbers)

    abundance_norm = abundances.sum()
    if abs(abundance_norm - 1) > 1e-12:
        logger.warn('Given abundances don\'t add up to 1 (value = %g) - normalizing', abundance_norm)
        abundances /= abundance_norm

    return abundances




def parse_spectral_bin(spectral_bin_boundary_1, spectral_bin_boundary_2):
    spectral_bin_boundary_1 = parse_quantity(spectral_bin_boundary_1).to('Angstrom', u.spectral())
    spectral_bin_boundary_2 = parse_quantity(spectral_bin_boundary_2).to('Angstrom', u.spectral())

    spectrum_start_wavelength = min(spectral_bin_boundary_1, spectral_bin_boundary_2)
    spectrum_end_wavelength = max(spectral_bin_boundary_1, spectral_bin_boundary_2)

    return spectrum_start_wavelength, spectrum_end_wavelength



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
            raise TARDISConfigurationError(
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
        raise TARDISConfigurationError('In abundance file section only types %s are allowed (supplied %s) ' %
                                (density_file_parser.keys(), density_file_dict['type']))

    return parser(density_file_dict, time_explosion)


def parse_velocity_section(velocity_dict, no_of_shells):
    velocity_parser = {}

    def parse_linear_velocity(velocity_dict, no_of_shells):
        v_inner = parse_quantity(velocity_dict['v_inner']).to('cm/s').value
        v_outer = parse_quantity(velocity_dict['v_outer']).to('cm/s').value
        velocities = np.linspace(v_inner, v_outer, no_of_shells + 1) * u.Unit('cm/s')
        return velocities[:-1], velocities[1:]

    velocity_parser['linear'] = parse_linear_velocity

    try:
        parser = velocity_parser[velocity_dict['type']]
    except KeyError:
        raise TARDISConfigurationError('In velocity section only types %s are allowed (supplied %s) ' %
                                (velocity_parser.keys(), velocity_dict['type']))
    return parser(velocity_dict, no_of_shells)


def parse_density_section(density_dict, no_of_shells, v_inner, v_outer, time_explosion):
    density_parser = {}


    #Parse density uniform
    def parse_uniform(density_dict, no_of_shells, v_inner, v_outer, time_explosion):

        return parse_quantity(density_dict['value']).to('g cm^-3') * np.ones(no_of_shells)

    density_parser['uniform'] = parse_uniform

    #Parse density branch85 w7
    def parse_branch85(density_dict, no_of_shells, v_inner, v_outer, time_explosion):

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

    def parse_exponential(density_dict, no_of_shells, v_inner, v_outer, time_explosion):
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
        raise TARDISConfigurationError('In density section only types %s are allowed (supplied %s) ' %
                                (density_parser.keys(), density_dict['type']))
    return parser(density_dict, no_of_shells, v_inner, v_outer, time_explosion)


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
        raise TARDISConfigurationError('In abundance file section only types %s are allowed (supplied %s) ' %
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

    def __dir__(self):
        return self.config_dict.keys()
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
    def from_yaml(cls, fname):
        try:
            yaml_dict = yaml.load(file(fname))
        except IOError as e:
            logger.critical('No config file named: %s', fname)
            raise e
        if yaml_dict['config_type'] not in ['simple1d']:
            raise TARDISConfigurationError('Only config_type=simple1d allowed at the moment.')

        return cls.from_config_dict(yaml_dict)

    @classmethod
    def from_config_dict(cls, raw_dict, atom_data=None):
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

        #First let's see if we can find an atom_db anywhere:

        if 'atom_data' in raw_dict.keys():
            atom_data_fname = raw_dict['atom_data']
        else:
            raise TARDISConfigurationError('No atom_data key found in config or command line')

        config_dict['atom_data_fname'] = atom_data_fname

        if atom_data is None:
            logger.info('Reading Atomic Data from %s', atom_data_fname)
            atom_data = atomic.AtomData.from_hdf5(atom_data_fname)
        else:
            atom_data = atom_data



        #Parsing supernova dictionary
        config_dict['supernova'] = parse_supernova_section(raw_dict['supernova'])

        #Trying to figure out the structure (number of shells)
        structure_section = raw_dict['model'].pop('structure')
        structure_config_dict = {}
        #first let's try to see if there's a file keyword
        if 'file' in structure_section:
            density_file_section = structure_section.pop('file')
            v_inner, v_outer, mean_densities, min_shell, max_shell = parse_density_file_section(density_file_section,
                                                                                                config_dict['supernova']
                                                                                                ['time_explosion'])

            no_of_shells = len(v_inner)
            if structure_section != {}:
                logger.warn(
                    'Accepted file for structure (density/velocity) structure ignoring all other arguments: \n%s\n',
                    pprint.pformat(structure_section, indent=4))
        else:
            #requiring all keys: no_of_shells, velocity, density
            if not all([item in structure_section.keys() for item in ('no_of_shells', 'velocity', 'density')]):
                raise TARDISConfigurationError(
                    'If file-section is not given to structure-section, one needs to provide all: no_of_shells, velocity, density')

            no_of_shells = structure_section['no_of_shells']

            v_inner, v_outer = parse_velocity_section(structure_section['velocity'], no_of_shells)
            mean_densities = parse_density_section(structure_section['density'], no_of_shells, v_inner, v_outer,
                                                   config_dict['supernova']['time_explosion'])

        r_inner = config_dict['supernova']['time_explosion'] * v_inner
        r_outer = config_dict['supernova']['time_explosion'] * v_outer
        r_middle = 0.5 * (r_inner + r_outer)


        structure_config_dict['v_inner'] = v_inner
        structure_config_dict['v_outer'] = v_outer
        structure_config_dict['mean_densities'] = mean_densities
        structure_config_dict['no_of_shells'] = no_of_shells
        structure_config_dict['r_inner'] = r_inner
        structure_config_dict['r_outer'] = r_outer
        structure_config_dict['r_middle'] = r_middle
        structure_config_dict['volumes'] = (4. / 3) * np.pi * (r_outer ** 3 - r_inner ** 3)




        config_dict['structure'] = structure_config_dict
        #Now that the structure section is parsed we move on to the abundances

        abundances_section = raw_dict['model']['abundances'].copy()
        abudances_config_dict = {}
        #TODO: columns are now until Z=120

        abundances = pd.DataFrame(columns=np.arange(no_of_shells),
                                  index=pd.Index(np.arange(1, 120), name='atomic_number'), dtype=np.float64)

        if 'file' in abundances_section.keys():
            abundance_file_dict = abundances_section.pop('file')
            parse_abundance_file_section(abundance_file_dict, abundances, min_shell, max_shell)


        for element_symbol_string in abundances_section:

            z = element_symbol2atomic_number(element_symbol_string, atom_data)

            abundances.ix[z] = float(abundances_section[element_symbol_string])

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

        if plasma_section['type'] not in ('nebular', 'lte'):
            raise TARDISConfigurationError('plasma_type only allowed to be "nebular" or "lte"')
        plasma_config_dict['type'] = plasma_section['type']

        if plasma_section['radiative_rates_type'] not in ('nebular', 'lte', 'detailed'):
            raise TARDISConfigurationError('radiative_rates_types must be either "nebular", "lte", or "detailed"')
        plasma_config_dict['radiative_rates_type'] = plasma_section['radiative_rates_type']

        if plasma_section['line_interaction_type'] not in ('scatter', 'downbranch', 'macroatom'):
            raise TARDISConfigurationError('radiative_rates_types must be either "scatter", "downbranch", or "macroatom"')
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
                    nlte_species.append(species_string_to_tuple(species_string, atom_data))

                nlte_config_dict['species'] = nlte_species
                nlte_config_dict['species_string'] = nlte_species_list
                nlte_config_dict.update(nlte_section)

                if 'coronal_approximation' not in nlte_section:
                    logger.debug('NLTE "coronal_approximation" not specified in NLTE section - defaulting to False')
                    nlte_config_dict['coronal_approximation'] = False

                if 'coronal_approximation' not in nlte_section:
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
        if 'last_no_of_packets' not in montecarlo_section:
            montecarlo_section['last_no_of_packets'] = None

        if 'no_of_virtual_packets' not in montecarlo_section:
            montecarlo_section['no_of_virtual_packets'] = 0

        montecarlo_config_dict.update(montecarlo_section)

        disable_electron_scattering = plasma_section['disable_electron_scattering']

        if disable_electron_scattering is False:
            logger.info("Electron scattering switched on")
            montecarlo_config_dict['sigma_thomson'] =6.652486e-25 / (u.cm**2)
        else:
            logger.warn('Disabling electron scattering - this is not physical')
            montecarlo_config_dict['sigma_thomson'] = 1e-200 / (u.cm**2)




        if 'black_body_sampling' in montecarlo_section:
            black_body_sampling_section = montecarlo_section.pop('black_body_sampling')
            sampling_start, sampling_end = parse_spectral_bin(black_body_sampling_section['start'],
                                                                                black_body_sampling_section['end'])
            montecarlo_config_dict['black_body_sampling']['start'] = sampling_start
            montecarlo_config_dict['black_body_sampling']['end'] = sampling_end
            montecarlo_config_dict['black_body_sampling']['samples'] = int(black_body_sampling_section['samples'])
        else:
            logger.warn('No "black_body_sampling" section in config file - using defaults of '
                        '50 - 200000 Angstrom (1e6 samples)')
            montecarlo_config_dict['black_body_sampling'] = {}
            montecarlo_config_dict['black_body_sampling']['start'] = 50 * u.angstrom
            montecarlo_config_dict['black_body_sampling']['end'] = 200000 * u.angstrom
            montecarlo_config_dict['black_body_sampling']['samples'] = int(1e6)

        config_dict['montecarlo'] = montecarlo_config_dict
        ##### End of MonteCarlo section






        ##### spectrum section ######
        spectrum_section = raw_dict.pop('spectrum')
        spectrum_config_dict = {}
        spectrum_start, spectrum_end = parse_spectral_bin(spectrum_section['start'], spectrum_section['end'])
        spectrum_bins = int(spectrum_section['bins'])
        spectrum_config_dict['start'] = spectrum_start
        spectrum_config_dict['end'] = spectrum_end
        spectrum_config_dict['bins'] = spectrum_bins

        spectrum_frequency = np.linspace(spectrum_end.to('Hz', u.spectral()).value,
                                                         spectrum_start.to('Hz', u.spectral()).value,
                                                         num=spectrum_bins+1)
        spectrum_config_dict['frequency'] = spectrum_frequency * u.Hz




        config_dict['spectrum'] = spectrum_config_dict


        return cls(config_dict, atom_data)


    def __init__(self, config_dict, atom_data):
        super(TARDISConfiguration, self).__init__(config_dict)
        self.atom_data = atom_data
        selected_atomic_numbers = self.abundances.index
        self.number_densities = (self.abundances * self.structure.mean_densities.to('g/cm^3').value)
        self.number_densities = self.number_densities.div(self.atom_data.atom_data.mass.ix[selected_atomic_numbers],
                                                          axis=0)








