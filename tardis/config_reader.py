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

logger = logging.getLogger(__name__)

data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))


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
    densities *= (time_explosion / time_0) ** -3

    return densities[1:]


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
    densities = rho_0 * (velocity_0 / velocities) ** exponent
    return  densities[1:]


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
    def from_file(cls, fname, args=None):
        config_parse_object = ConfigParser()
        config_parse_object.read(fname)
        general_dict = dict(config_parse_object.items('general'))
        abundance_dict = dict(config_parse_object.items('abundances'))

        config_object = cls()
        config_object.parse_general_section(general_dict)
        config_object.parse_abundance_section(abundance_dict)
        return config_object


    def __init__(self):
        self.velocities = None
        self.densities = None
        self.abundances = None
        self.initial_t_rad = 10000
        self.t_inner = 10000.0  # TODO calculate
        self.ws = None
        self.no_of_shells = None

        self._luminosity_outer = None
        self.time_of_simulation = None


    def parse_general_section(self, config_dict):
        model_type = config_dict.pop('model_type')

        if model_type != 'radial1d':
            raise ValueError("Only supporting 'radial1d' at the moment")

        # reading time since explosion
        time_explosion_value, time_explosion_unit = config_dict.pop('time_explosion').split()
        self.time_explosion = units.Quantity(float(time_explosion_value), time_explosion_unit).to('s').value

        # Reading luminosity, special unit log_l_sun is luminosity given in log10
        # of solar units
        luminosity_value, luminosity_unit = config_dict.pop('luminosity').split()
        if luminosity_unit == 'log_lsun':
            self.luminosity_outer = 10 ** (
                float(luminosity_value) + np.log10(constants.cgs.L_sun.value))
        else:
            self.luminosity_outer = units.Quantity(
                float(luminosity_value), luminosity_unit).to('erg/s').value

        # reading number of shells
        no_of_shells = int(config_dict.pop('zones'))

        self.no_of_shells = no_of_shells

        # reading velocities
        # set of velocities currently supported are v_inner, v_outer and
        # v_sampling linear

        v_inner_value, v_inner_unit = config_dict.pop('v_inner').split()
        v_inner = units.Quantity(
            float(v_inner_value), v_inner_unit).to('cm/s').value

        v_outer_value, v_outer_unit = config_dict.pop('v_outer').split()
        v_outer = units.Quantity(
            float(v_outer_value), v_outer_unit).to('cm/s').value

        v_sampling = config_dict.pop('v_sampling')

        self.set_velocities(
            v_inner=v_inner, v_outer=v_outer, v_sampling=v_sampling)

        density_set = config_dict.pop('density')

        if density_set == 'w7_branch85':
            self.densities = calculate_w7_branch85_densities(
                self.velocities,
                self.time_explosion)
        elif density_set == 'exponential':
            #TODO:Add here the function call which generates the exponential density profile. The easy way from tonight don't  work as expected!!
            if not (('exponential_n_factor' in config_dict) and ('exponential_rho0' in config_dict)):
                raise ValueError(
                    'If density=exponential is set the exponential_n_factor(float) and exponential_rho_0 have to be specified.')

            self.exponential_n_factor = float(config_dict.pop('exponential_n_factor'))
            self.exponential_rho_0 = float(config_dict.pop('exponential_rho0'))

            self.densities = calculate_exponential_densities(self.velocities, v_inner,
                self.exponential_rho_0, self.exponential_n_factor)


        else:
            try:
                density = float(density_set)
            except ValueError:
                raise ValueError(
                    'Currently only density = w7_branch85 or density = exponential,'
                    ' or specifying a uniform density (with a number) are supported')
            self.densities = np.ones(self.no_of_shells) * density


        # reading plasma type
        self.plasma_type = config_dict.pop('plasma_type')

        # reading initial t_rad
        if 'initial_t_rad' in config_dict:
            initial_t_rad_value, initial_t_rad_unit = config_dict.pop('initial_t_rad').split()
            self.initial_t_rad = units.Quantity(float(initial_t_rad_value), initial_t_rad_unit).to('K').value
            logger.info('Selected %g K as initial temperature', self.initial_t_rad)
        else:
            logger.warn('No initial shell temperature specified (initial_t_rad) - using default 10000 K')

        # reading line interaction type
        self.line_interaction_type = config_dict.pop('line_interaction_type')

        atom_data_file = config_dict.pop('atom_data_file', None)

        if atom_data_file is None:
            raise ValueError("Please specify a filename with the keyword 'atom_data_file'")

        if self.plasma_type.lower == 'lte':
            self.atom_data = atomic.AtomData.from_hdf5(atom_data_file)
        else:
            self.atom_data = atomic.AtomData.from_hdf5(atom_data_file, use_macro_atom=True,
                use_zeta_data=True)


        # reading number of packets and iterations
        if 'iterations' in config_dict and 'no_of_packets' in config_dict:
            self.iterations = int(float(config_dict.pop('iterations')))
            self.no_of_packets = int(float(config_dict.pop('no_of_packets')))
        else:
            raise ValueError("'no_of_packets' and 'iterations' needs to be set in the configuration")

        last_no_of_packets = config_dict.pop('last_no_of_packets', None)
        if last_no_of_packets is not None:
            self.last_no_of_packets = int(float(last_no_of_packets))
            logger.info('Last iteration will have %g packets', self.last_no_of_packets)

        no_of_virtual_packets = config_dict.pop('no_of_virtual_packets', None)

        if no_of_virtual_packets is not None:
            self.no_of_virtual_packets = int(float(no_of_virtual_packets))
            logger.info('Activating Virtual packets for last iteration (%g)', self.no_of_virtual_packets)

        spectrum_start_value, spectrum_end_unit = config_dict.pop(
            'spectrum_start').split()
        spectrum_start = units.Quantity(float(spectrum_start_value), spectrum_end_unit).to('angstrom',
            units.spectral()).value

        spectrum_end_value, spectrum_end_unit = config_dict.pop('spectrum_end').split()
        spectrum_end = units.Quantity(float(spectrum_end_value), spectrum_end_unit).to('angstrom',
            units.spectral()).value

        self.spectrum_bins = int(float(config_dict.pop('spectrum_bins')))

        if spectrum_end > spectrum_start:
            logger.debug('Converted spectrum start/end to angstrom %.4g %.4g', spectrum_start, spectrum_end)
            self.spectrum_start = spectrum_start
            self.spectrum_end = spectrum_end

        else:
            logger.warn('Spectrum Start > Spectrum End in wavelength space - flipped them')

            logger.debug('Converted spectrum start/end to angstrom %.4g %.4g', spectrum_end, spectrum_start)

            self.spectrum_start = spectrum_end
            self.spectrum_end = spectrum_start

        self.spectrum_start_nu = units.Quantity(self.spectrum_end, 'angstrom').to('Hz', units.spectral())
        self.spectrum_end_nu = units.Quantity(self.spectrum_start, 'angstrom').to('Hz', units.spectral())

        if config_dict != {}:
            logger.warn('Not all config options parsed - ignored %s' % config_dict)

    def parse_abundance_section(self, abundance_dict):
        abundances = {}
        abundance_set = abundance_dict.pop('abundance_set', None)
        if abundance_set == 'lucy99':
            abundances = read_lucy99_abundances()
        elif abundance_set is not None:
            raise ValueError("Currently only 'lucy99' abundance_set supported")

        self.nlte_species = []
        species_pattern = re.compile('\s*([a-zA-Z]*)(\d*)\s*')
        if 'nlte_species' in abundance_dict:
            for species_symbol in abundance_dict.pop('nlte_species').split(','):
                species_match = species_pattern.match(species_symbol)
                if species_match is None:
                    raise ValueError(
                        "'nlte_species' element %s could not be matched to a valid species notation (e.g. Ca2)")
                species_element, species_ion = species_match.groups()
                species_element = reformat_element_symbol(species_element)
                if species_element not in self.atom_data.symbol2atomic_number:
                    raise ValueError("Element provided in NLTE species %s unknown" % species_element)
                self.nlte_species.append((self.atom_data.symbol2atomic_number[species_element], int(species_ion) - 1))

        for element in abundance_dict:
            element_symbol = reformat_element_symbol(element)
            if element_symbol not in self.atom_data.symbol2atomic_number:
                raise ValueError('Element %s provided in config unknown' % element_symbol)
            if element_symbol in abundances:
                logger.debug('Element %s already in abundances - overwriting %g with %g', (element_symbol,
                                                                                           abundances[element_symbol],
                                                                                           abundance_dict[element]))
            abundances[element_symbol] = float(abundance_dict[element])

        self.set_abundances(abundances)


    @property
    def line_interaction_type(self):
        return self._line_interaction_type

    @line_interaction_type.setter
    def line_interaction_type(self, value):
        if value not in ('scatter', 'downbranch', 'macroatom'):
            raise ValueError('line_interaction_type can only be "scatter", "downbranch", or "macroatom"')
        self._line_interaction_type = value


    def set_velocities(self, velocities=None, v_inner=None, v_outer=None, v_sampling='linear'):
        """
        Setting the velocities

        :param velocities:
        :param v_inner:
        :param v_outer:
        :param v_sampling:


        """
        if self.no_of_shells is None:
            raise ValueError('Can not set abundances before number of shells have been set')

        if v_sampling == 'linear':
            self.velocities = np.linspace(
                v_inner, v_outer, self.no_of_shells + 1)
        else:
            raise ValueError('Currently only v_sampling = linear is possible')

    def set_abundances(self, abundances):
        """
        Setting the abundances

        abundances: `dict` or `list`
            if a dict is given the assumed mode is uniform, if a list is given it must contain only lists


        """
        if self.no_of_shells is None:
            raise ValueError('Can not set abundances before number of shells have been set')

        if isinstance(abundances, dict):
            self.abundances = [abundances] * self.no_of_shells

    @property
    def luminosity_outer(self):
        return self._luminosity_outer

    @luminosity_outer.setter
    def luminosity_outer(self, value):
        self._luminosity_outer = value

    def set_densities(self, densities):
        """

        :param densities:
        :return:
        """

        self.densities = densities

    def final_preparation(self):
        """
        Does the final preparation for the configuration object

        Generates the atoms needed in the simulation

        :return:
        """
        self.selected_atoms = set()
        for shell_abundance in self.abundances:
            self.selected_atoms.update(shell_abundance.keys())

            # self.required_atomic_number =
            # set([atom_data.symbol2atomic_number[item] for item in
            # required_atoms])


