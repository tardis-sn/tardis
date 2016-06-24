# This module contains the model class

import logging
import os

import numpy as np
import pandas as pd
from astropy import constants, units as u
from astropy.utils.decorators import deprecated

from tardis.io.util import to_hdf
from util import intensity_black_body

from tardis.plasma.standard_plasmas import LegacyPlasmaArray

logger = logging.getLogger(__name__)

c = constants.c.cgs.value
h = constants.h.cgs.value
kb = constants.k_B.cgs.value


class Radial1DModel(object):
    """
        Class to hold the states of the individual shells (the state of the plasma (as a `~plasma.BasePlasma`-object or one of its subclasses),
        , the plasma parameters (e.g. temperature, dilution factor), the dimensions of the shell).


        Parameters
        ----------

        tardis_configuration : `tardis.config_reader.Configuration`

        velocities : `np.ndarray`
            an array with n+1 (for n shells) velocities (in cm/s) for each of the boundaries (velocities[0] describing
            the inner boundary and velocities[-1] the outer boundary

        densities : `np.ndarray`
            an array with n densities - being the density mid-shell (assumed for the whole shell)

        abundances : `list` or `dict`
            a dictionary for uniform abundances throughout all shells, e.g. dict(Fe=0.5, Si=0.5)
            For a different abundance for each shell list of abundance dictionaries.


        time_explosion : `float`
            time since explosion in seconds

        atom_data : `~tardis.atom_data.AtomData` class or subclass
            Containing the atom data needed for the plasma calculations

        ws : `None` or `list`-like
            ws can only be specified for plasma_type 'nebular'. If `None` is specified at first initialization the class
            calculates an initial geometric dilution factor. When giving a list positive values will be accepted, whereas
            negative values trigger the usage of the geometric calculation

        plasma_type : `str`
            plasma type currently supports 'lte' (using `tardis.plasma.LTEPlasma`)
            or 'nebular' (using `tardis.plasma.NebularPlasma`)

        initial_t_rad : `float`-like or `list`-like
            initial radiative temperature for each shell, if a scalar is specified it initializes with a uniform
            temperature for all shells




    """

    @classmethod
    def from_h5(cls, buffer_or_fname):
        raise NotImplementedError("This is currently not implemented")


    def __init__(self, tardis_config):
        #final preparation for configuration object
        self.tardis_config = tardis_config

        self.atom_data = tardis_config.atom_data
        selected_atomic_numbers = self.tardis_config.abundances.index

        self.atom_data.prepare_atom_data(
            selected_atomic_numbers,
            line_interaction_type=tardis_config.plasma.line_interaction_type,
            nlte_species=tardis_config.plasma.nlte.species)

        if tardis_config.plasma.ionization == 'nebular':
            if not self.atom_data.has_zeta_data:
                raise ValueError("Requiring Recombination coefficients Zeta "
                                 "for 'nebular' plasma ionization")

        self.t_inner = tardis_config.plasma.t_inner

        self.v_inner = tardis_config.structure.v_inner
        self.v_outer = tardis_config.structure.v_outer
        self.v_middle = 0.5 * (self.v_inner + self.v_outer)

        self.ws = self.calculate_geometric_w(
            tardis_config.structure.r_middle,
            tardis_config.structure.r_inner[0])

        if tardis_config.plasma.t_rads is None:
            self.t_rads = self._init_t_rad(
                self.t_inner, self.v_inner[0], self.v_middle)
        else:
            self.t_rads = tardis_config.plasma.t_rads

        heating_rate_data_file = getattr(
            tardis_config.plasma, 'heating_rate_data_file', None)

        self.plasma = LegacyPlasmaArray(
            tardis_config.number_densities, tardis_config.atom_data,
            tardis_config.supernova.time_explosion.to('s').value,
            nlte_config=tardis_config.plasma.nlte,
            delta_treatment=tardis_config.plasma.get('delta_treatment', None),
            ionization_mode=tardis_config.plasma.ionization,
            excitation_mode=tardis_config.plasma.excitation,
            line_interaction_type=tardis_config.plasma.line_interaction_type,
            link_t_rad_t_electron=0.9,
            helium_treatment=tardis_config.plasma.helium_treatment,
            heating_rate_data_file=heating_rate_data_file,
            v_inner=self.v_inner,
            v_outer=self.v_outer)

        distance = tardis_config.supernova.get('distance', None)
        self.spectrum = TARDISSpectrum(
            tardis_config.spectrum.frequency, distance)
        self.spectrum_virtual = TARDISSpectrum(
            tardis_config.spectrum.frequency, distance)
        self.spectrum_reabsorbed = TARDISSpectrum(
            tardis_config.spectrum.frequency, distance)

        self.calculate_j_blues(init_detailed_j_blues=True)
        self.update_plasmas(initialize_nlte=True)


    @property
    @deprecated('v1.5',
                'plasma_array has been renamed to plasma and will be removed in the future. Please use model.plasma instead.')
    def plasma_array(self):
        return self.plasma

    @plasma_array.setter
    def plasma_array(self, value):
        self.plasma = value

    @property
    def line_interaction_type(self):
        return self._line_interaction_type

    @line_interaction_type.setter
    def line_interaction_type(self, value):
        if value in ['scatter', 'downbranch', 'macroatom']:
            self._line_interaction_type = value
            self.tardis_config.plasma.line_interaction_type = value
            #final preparation for atom_data object - currently building data
            self.atom_data.prepare_atom_data(
                self.tardis_config.number_densities.columns,
                line_interaction_type=self.line_interaction_type,
                max_ion_number=None,
                nlte_species=self.tardis_config.plasma.nlte.species)
        else:
            raise ValueError('line_interaction_type can only be '
                             '"scatter", "downbranch", or "macroatom"')



    @property
    def t_inner(self):
        return self._t_inner

    @t_inner.setter
    def t_inner(self, value):
        self._t_inner = value
        self.luminosity_inner = (
            4 * np.pi * constants.sigma_sb.cgs *
            self.tardis_config.structure.r_inner[0] ** 2
            * self.t_inner ** 4).to('erg/s')

        self.time_of_simulation = (1.0 * u.erg / self.luminosity_inner)
        self.j_blues_norm_factor = (
            constants.c.cgs * self.tardis_config.supernova.time_explosion /
            (4 * np.pi * self.time_of_simulation *
             self.tardis_config.structure.volumes))


    @staticmethod
    def calculate_geometric_w(r, r_inner):
        return 0.5 * (1 - np.sqrt(1 - (r_inner ** 2 / r ** 2).to(1).value))

    @staticmethod
    def _init_t_rad(t_inner, v_boundary, v_middle):
        lambda_wien_inner = constants.b_wien / t_inner
        return constants.b_wien / (
            lambda_wien_inner * (1 + (v_middle - v_boundary) / constants.c))


    def calculate_j_blues(self, init_detailed_j_blues=False):
        nus = self.atom_data.lines.nu.values
        radiative_rates_type = self.tardis_config.plasma.radiative_rates_type
        w_epsilon = self.tardis_config.plasma.w_epsilon

        if radiative_rates_type == 'blackbody':
            logger.info('Calculating J_blues for radiative_rates_type=lte')
            j_blues = intensity_black_body(nus[np.newaxis].T, self.t_rads.value)
            self.j_blues = pd.DataFrame(
                j_blues, index=self.atom_data.lines.index,
                columns=np.arange(len(self.t_rads)))

        elif radiative_rates_type == 'dilute-blackbody' or init_detailed_j_blues:
            logger.info('Calculating J_blues for radiative_rates_type=dilute-blackbody')
            j_blues = self.ws * intensity_black_body(nus[np.newaxis].T, self.t_rads.value)
            self.j_blues = pd.DataFrame(
                j_blues, index=self.atom_data.lines.index,
                columns=np.arange(len(self.t_rads)))

        elif radiative_rates_type == 'detailed':
            logger.info('Calculating J_blues for radiate_rates_type=detailed')

            self.j_blues = pd.DataFrame(
                self.j_blue_estimators *
                    self.j_blues_norm_factor.value,
                    index=self.atom_data.lines.index,
                    columns=np.arange(len(self.t_rads)))

            for i in xrange(self.tardis_config.structure.no_of_shells):
                zero_j_blues = self.j_blues[i] == 0.0
                self.j_blues[i][zero_j_blues] = (
                    w_epsilon * intensity_black_body(
                        self.atom_data.lines.nu[zero_j_blues].values,
                        self.t_rads.value[i]))

        else:
            raise ValueError('radiative_rates_type type unknown - %s', radiative_rates_type)

    def update_plasmas(self, initialize_nlte=False):

        self.plasma.update_radiationfield(
            self.t_rads.value, self.ws, self.j_blues,
            self.tardis_config.plasma.nlte, initialize_nlte=initialize_nlte,
            n_e_convergence_threshold=0.05)

        if self.tardis_config.plasma.line_interaction_type in ('downbranch',
                                                               'macroatom'):
            self.transition_probabilities = (
                self.plasma.transition_probabilities)

    def save_spectra(self, fname):
        self.spectrum.to_ascii(fname)
        self.spectrum_virtual.to_ascii('virtual_' + fname)

    def to_hdf(self, path_or_buf, path='', plasma_properties=None):
        """
        Store the model to an HDF structure.

        Parameters
        ----------
        path_or_buf
            Path or buffer to the HDF store
        path : str
            Path inside the HDF store to store the model
        plasma_properties
            `None` or a `PlasmaPropertyCollection` which will
            be passed as the collection argument to the
            plasma.to_hdf method.

        Returns
        -------
        None

        """
        model_path = os.path.join(path, 'model')
        properties = ['t_inner', 'ws', 't_rads', 'v_inner', 'v_outer']
        to_hdf(path_or_buf, model_path, {name: getattr(self, name) for name
                                         in properties})

        self.plasma.to_hdf(path_or_buf, model_path, plasma_properties)

        metadata = pd.Series({'atom_data_uuid': self.atom_data.uuid1})
        metadata.to_hdf(path_or_buf,
                                 os.path.join(model_path, 'metadata'))



class TARDISSpectrum(object):
    """
    TARDIS Spectrum object
    """

    def __init__(self, frequency, distance=None):
        self._frequency = frequency
        self.wavelength = self.frequency.to('angstrom', u.spectral())

        self.distance = distance



        self.delta_frequency = frequency[1] - frequency[0]

        self._flux_nu = np.zeros_like(frequency.value) * u.Unit('erg / (s Hz cm^2)')
        self._flux_lambda = np.zeros_like(frequency.value) * u.Unit('erg / (s Angstrom cm^2)')

        self.luminosity_density_nu = np.zeros_like(self.frequency) * u.Unit('erg / (s Hz)')
        self.luminosity_density_lambda = np.zeros_like(self.frequency) * u.Unit('erg / (s Angstrom)')

    @property
    def frequency(self):
        return self._frequency[:-1]


    @property
    def flux_nu(self):
        if self.distance is None:
            raise AttributeError('supernova distance not supplied - flux calculation impossible')
        else:
            return self._flux_nu

    @property
    def flux_lambda(self):
        if self.distance is None:
            raise AttributeError('supernova distance not supplied - flux calculation impossible')
        return self._flux_lambda

    def update_luminosity(self, spectrum_luminosity):
        self.luminosity_density_nu = (spectrum_luminosity / self.delta_frequency).to('erg / (s Hz)')
        self.luminosity_density_lambda = self.f_nu_to_f_lambda(self.luminosity_density_nu.value) \
                                         * u.Unit('erg / (s Angstrom)')

        if self.distance is not None:
            self._flux_nu = (self.luminosity_density_nu / (4 * np.pi * self.distance.to('cm')**2))

            self._flux_lambda = self.f_nu_to_f_lambda(self.flux_nu.value) * u.Unit('erg / (s Angstrom cm^2)')

    def f_nu_to_f_lambda(self, f_nu):
        return f_nu * self.frequency.value**2 / constants.c.cgs.value / 1e8


    def plot(self, ax, mode='wavelength'):
        if mode == 'wavelength':
            ax.plot(self.wavelength.value, self.flux_lambda.value)
            ax.set_xlabel('Wavelength [%s]' % self.wavelength.unit._repr_latex_())
            ax.set_ylabel('Flux [%s]' % self.flux_lambda.unit._repr_latex_())

    def to_ascii(self, fname, mode='luminosity_density'):
        if mode == 'luminosity_density':
            np.savetxt(fname, zip(self.wavelength.value, self.luminosity_density_lambda.value))
        elif mode == 'flux':
            np.savetxt(fname, zip(self.wavelength.value, self.flux_lambda.value))
        else:
            raise NotImplementedError('only mode "luminosity_density" and "flux" are implemented')
