# Module to read the rather complex config data
# Currently the configuration file is documented in
# tardis/data/example_configuration.ini

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


def calculate_exponentail_densities(velocities, velocity_0, rho_0, exponent):
    """

    :param velocities: Array like list of velocities
    :param velocity_0: Velocity at the inner boundary
    :param rho_0: density at velocity_0
    :param exponent: exponent used in the powerlaw
    :return: Array like density structure

    This function computes a descret exponential density profile.
    math: \rho = \rho_0 \times \left( \frac{v_0}{v} \right)^n
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


def read_config(fname, atomic_data=None):
    config_object = ConfigParser()
    config_object.read(fname)
    tardis_configuration = TardisConfiguration()
    general_dict = dict(config_object.items('general'))
    parse_general_section(general_dict, tardis_configuration)
    abundance_dict = dict(config_object.items('abundances'))
    tardis_configuration.set_abundances(
        parse_abundance_section(abundance_dict, atomic_data))
    return tardis_configuration


class TardisConfiguration(object):
    """
    Tardis configuration class
    """

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
        self.exponential_n_factor = 10
        self.exponential_rho_0 = 10e5
        self.simulation_type = 'single_run'

    @property
    def number_of_packets(self):
        """

        returns the right number of packets with single run and will keep track of iterations with
        multirun

        """
        if self.simulation_type == 'single_run':
            return self.single_run_packets
        else:
            raise ValueError('Currently only single_run supported')

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


def parse_abundance_section(abundance_dict, atomic_data=None ):
    if atomic_data is None: # fallback in case no atomic dataset is given!
        atomic_dict = {"H": "1"}
        logger.warn('Using fallback because no atomic dataset is given')
    else:
        atomic_dict = dict(atomic_data.symbol2atomic_number)

    abundance_set = abundance_dict.get('abundance_set', None)

    if abundance_set == 'lucy99':
        abundances = read_lucy99_abundances()
    else:
        #Added by Michi
        abundance_set = dict()
        print(atomic_dict.values())

        for current_name in atomic_dict.keys():
            current_cont = abundance_dict.get(current_name.lower(), 0.)
            #print(current_cont)
            print(current_name)
            if current_cont != 0:
                abundance_set[current_name] = float(current_cont)

        print(abundance_set.values())

        abundances_sum = sum(abundance_set.values())
        if abundances_sum < 1.:
            abundance_set['H'] = 1 - abundances_sum
        elif  abundances_sum > 1:
            for current_name in abundance_set:
                abundance_set[current_name] *= 1. / abundances_sum


        #raise ValueError('Currently only abundance_set=lucy99 supported')

        abundances = abundance_set
    return abundances


def parse_general_section(config_dict, general_config):
    model_type = config_dict.pop('model_type')

    if model_type != 'radial1d':
        raise ValueError("Only supporting 'radial1d' at the moment")

    # reading time since explosion
    time_explosion_value, time_explosion_unit = config_dict.pop(
        'time_explosion').split()
    general_config.time_explosion = units.Quantity(
        float(time_explosion_value), time_explosion_unit).to('s').value

    # Reading luminosity, special unit log_l_sun is luminosity given in log10
    # of solar units
    luminosity_value, luminosity_unit = config_dict.pop('luminosity').split()
    if luminosity_unit == 'log_lsun':
        general_config.luminosity_outer = 10 ** (
            float(luminosity_value) + np.log10(constants.cgs.L_sun.value))
    else:
        general_config.luminosity_outer = units.Quantity(
            float(luminosity_value), luminosity_unit).to('erg/s').value

    # reading number of shells
    no_of_shells = int(config_dict.pop('zones'))

    general_config.no_of_shells = no_of_shells

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

    general_config.set_velocities(
        v_inner=v_inner, v_outer=v_outer, v_sampling=v_sampling)

    density_set = config_dict.pop('density_set')

    if density_set == 'w7_branch85':
        general_config.densities = calculate_w7_branch85_densities(
            general_config.velocities,
            general_config.time_explosion)
    elif density_set == 'exponential':
        #TODO:Add here the function call which generates the exponential density profile. The easy way from tonight don't  work as expected!!
        if (('exponential_n_factor' in config_dict) & ('exponential_rho0' in config_dict)):
            try:
                general_config.exponential_n_factor = float(config_dict.pop('exponential_n_factor'))
                general_config.exponential_rho_0 = float(config_dict.pop('exponential_rho0'))
            except ValueError:
                logger.warn(
                    'If density_set=exponential is set the exponential_n_factor(float) and exponential_rho_0 have to be specified. Using the default density_set=10! ')

            general_config.densities = calculate_exponentail_densities(general_config.velocities, v_inner,
                general_config.exponential_rho_0, general_config.exponential_n_factor)
            print(general_config.densities)
    else:
        raise ValueError(
            'Curently only density_set = w7_branch85 or density_set = exponential are supported')

    # reading plasma type
    general_config.plasma_type = config_dict.pop('plasma_type')

    # reading initial t_rad
    if 'initial_t_rad' in config_dict:
        general_config.initial_t_rad = float(config_dict.pop('initial_t_rad'))
    else:
        logger.warn('No initial shell temperature specified (initial_t_rad) - using default 10000 K')

    # reading line interaction type
    general_config.line_interaction_type = config_dict.pop(
        'line_interaction_type')

    # reading number of packets and iterations
    if 'single_run_packets' in config_dict:
        if [item for item in ('spectrum_packets', 'calibration_packets') if item in config_dict]:
            raise ValueError('Please specify either "spectrum_packets"/"calibration_packets"/"iterations" or '
                             '"single_run_packets" in config file')

        general_config.simulation_type = 'single_run'
        general_config.single_run_packets = int(
            float(config_dict.pop('single_run_packets')))
    elif len([item for item in ('spectrum_packets', 'calibration_packets', 'iterations') if item in config_dict]) == 3:
        general_config.calibration_packets = int(
            float(config_dict.pop('calibration_packets')))
        general_config.spectrum_packets = int(
            float(config_dict.pop('spectrum_packets')))
        general_config.iterations = int(float(config_dict.pop('iterations')))
    else:
        raise ValueError('Please specify either "spectrum_packets"/"calibration_packets"/"iterations" or'
                         ' "single_run_packets" in config file')

    # TODO fix quantity spectral in astropy

    spectrum_start_value, spectrum_start_unit = config_dict.pop(
        'spectrum_start').split()
    general_config.spectrum_start_wave = units.Quantity(
        float(spectrum_start_value), spectrum_start_unit).value
    general_config.spectrum_end_nu = units.Unit('angstrom').to('Hz', general_config.spectrum_start_wave, units.spectral())
    spectrum_end_value, spectrum_end_unit = config_dict.pop(
        'spectrum_end').split()
    general_config.spectrum_end_wave = units.Quantity(
        float(spectrum_end_value), spectrum_end_unit).value
    general_config.spectrum_start_nu = units.Unit('angstrom').to('Hz', general_config.spectrum_end_wave, units.spectral())
 
    general_config.spectrum_bins = int(float (config_dict.pop('spectrum_bins')))


    if config_dict != {}:
        logger.warn('Not all config options parsed - ignored %s' % config_dict)

    return general_config
