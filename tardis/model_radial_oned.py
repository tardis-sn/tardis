# building of radial_oned_model

import numpy as np
import plasma, atomic, packet_source
import logging
import config_reader
import pandas as pd
from astropy import constants
from copy import deepcopy

logger = logging.getLogger(__name__)

c = constants.cgs.c.value
h = constants.cgs.h.value
kb = constants.cgs.k_B.value

w_estimator_constant = (c ** 2 / (2 * h)) * (15 / np.pi ** 4) * (h / kb) ** 4 / (4 * np.pi)


class Radial1DModel(object):
    """
        Class to hold the states of the individual shells (the state of the plasma (as a `~plasma.BasePlasma`-object or one of its subclasses),
        , the plasma parameters (e.g. temperature, dilution factor), the dimensions of the shell).


        Parameters
        ----------

        tardis_configuration : `tardis.config_reader.TardisConfiguration`

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


    def __init__(self, tardis_config):
        #final preparation for configuration object
        tardis_config.final_preparation()

        self.tardis_config = tardis_config

        self.atom_data = tardis_config.atom_data

        self.packet_src = packet_source.SimplePacketSource.from_wavelength(tardis_config.spectrum_start,
            tardis_config.spectrum_end)

        self.no_of_shells = tardis_config.no_of_shells

        logger.info('Assuming %d shells' % self.no_of_shells)


        #setting time_explosion
        self.time_explosion = tardis_config.time_explosion

        #initializing velocities and radii
        self.v_inner = tardis_config.velocities[:-1]
        self.v_outer = tardis_config.velocities[1:]

        self.r_inner = self.v_inner * self.time_explosion
        self.r_outer = self.v_outer * self.time_explosion
        self.r_middle = 0.5 * (self.r_inner + self.r_outer)

        self.volumes = (4. / 3) * np.pi * (self.r_outer ** 3 - self.r_inner ** 3)
        #initializing densities
        assert len(tardis_config.densities) == self.no_of_shells

        self.densities_middle = tardis_config.densities

        self.time_of_simulation = tardis_config.time_of_simulation

        self.t_inner = tardis_config.t_inner

        self.luminosity_outer = tardis_config.luminosity_outer

        self.no_of_packets = tardis_config.no_of_packets
        self.iterations = tardis_config.iterations
        self.create_packets()

        self.spec_virt_nu = np.linspace(configuration_object.spectrum_start_nu, configuration_object.spectrum_end_nu, configuration_object.spectrum_bins+1)
        
        self.spec_virt_flux_nu = np.zeros_like(self.spec_virt_nu)
        
        
    #Selecting plasma class
        self.plasma_type = configuration_object.plasma_type

        #Selecting plasma class
        self.plasma_type = tardis_config.plasma_type
        if self.plasma_type == 'lte':
            self.plasma_class = plasma.LTEPlasma
            if tardis_config.ws is not None:
                raise ValueError(
                    "the dilution factor W ('ws') can only be specified when selecting plasma_type='nebular'")

        elif self.plasma_type == 'nebular':
            self.plasma_class = plasma.NebularPlasma
            if not self.atom_data.has_zeta_data:
                raise ValueError("Requiring Recombination coefficients Zeta for 'nebular' plasma_type")
        else:
            raise ValueError("Currently this model only supports 'lte' or 'nebular'")




        #initializing abundances
        self.abundances = tardis_config.abundances
        self.number_densities = []
        for abundance, density in zip(self.abundances, self.densities_middle):
            self.number_densities.append(calculate_atom_number_densities(self.atom_data, abundance, density))

        self.selected_atomic_numbers = self.number_densities[0].index.values.astype(int)

        self.line_interaction_type = tardis_config.line_interaction_type


        #setting dilution factors
        if tardis_config.ws is None:
            self.ws = 0.5 * (1 - np.sqrt(1 - self.r_inner[0] ** 2 / self.r_middle ** 2))
        else:
            self.ws = np.array([(0.5 * (1 - np.sqrt(1 - self.r_inner[0] ** 2 / self.r_middle[i] ** 2))) if w < 0\
                                else w for i, w in enumerate(tardis_config.ws)])

        #initializing temperatures

        if np.isscalar(tardis_config.initial_t_rad):
            self.t_rads = [tardis_config.initial_t_rad] * self.no_of_shells
        else:
            assert len(tardis_config.initial_t_rad) == self.no_of_shells
            self.t_rads = np.array(tardis_config.initial_t_rad, dtype=np.float64)

        self.initialize_plasmas()


    @property
    def electron_density(self):
        return np.array([plasma.electron_density for plasma in self.plasmas])

    @property
    def line_interaction_type(self):
        return self._line_interaction_type

    @line_interaction_type.setter
    def line_interaction_type(self, value):
        if value == 'scatter':
            self.line_interaction_id = 0
        elif value == 'downbranch':
            self.line_interaction_id = 1
        elif value == 'macroatom':
            self.line_interaction_id = 2
        else:
            raise ValueError('line_interaction_type can only be "scatter", "downbranch", or "macroatom"')

        self._line_interaction_type = value
        #final preparation for atom_data object - currently building data
        self.atom_data.prepare_atom_data(self.selected_atomic_numbers,
            line_interaction_type=self.line_interaction_type, max_ion_number=None)

    @property
    def t_inner(self):
        return self._t_inner

    @t_inner.setter
    def t_inner(self, value):
        self._t_inner = value
        self.luminosity_inner = 4 * np.pi * constants.cgs.sigma_sb.value * self.r_inner[0] ** 2 * self._t_inner ** 4
        self.time_of_simulation = 1 / self.luminosity_inner


    def create_packets(self):
        #Energy emitted from the inner boundary
        self.emitted_inner_energy = 4 * np.pi * constants.cgs.sigma_sb.value * self.r_inner[0] ** 2 * (
            self.t_inner) ** 4
        self.packet_src.create_packets(self.no_of_packets, self.t_inner)

    def initialize_plasmas(self):
        self.plasmas = []
        self.tau_sobolevs = np.empty((self.no_of_shells, len(self.atom_data.lines)))
        self.line_list_nu = self.atom_data.lines['nu']

        if self.line_interaction_id in (1, 2):
            self.transition_probabilities = []

        if self.plasma_type == 'lte':
            for i, (current_abundances, current_t_rad) in\
            enumerate(zip(self.number_densities, self.t_rads)):
                current_plasma = self.plasma_class(current_abundances, self.atom_data, self.time_explosion,
                    nlte_species=self.tardis_config.nlte_species)
                current_plasma.update_radiationfield(current_t_rad)
                self.tau_sobolevs[i] = current_plasma.tau_sobolevs

                #TODO change this
                if self.line_interaction_id in (1, 2):
                    self.transition_probabilities.append(current_plasma.update_macro_atom().values)

                self.plasmas.append(current_plasma)


        elif self.plasma_type == 'nebular':
            for i, (current_abundances, current_t_rad, current_w) in\
            enumerate(zip(self.number_densities, self.t_rads, self.ws)):
                current_plasma = self.plasma_class(current_abundances, self.atom_data, self.time_explosion,
                    nlte_species=self.tardis_config.nlte_species)
                current_plasma.update_radiationfield(current_t_rad, current_w)
                self.tau_sobolevs[i] = current_plasma.tau_sobolevs

                if self.line_interaction_id in (1, 2):
                    self.transition_probabilities.append(current_plasma.update_macro_atom().values)

                self.plasmas.append(current_plasma)

        self.tau_sobolevs = np.array(self.tau_sobolevs, dtype=float)

        if self.line_interaction_id in (1, 2):
            self.transition_probabilities = np.array(self.transition_probabilities, dtype=np.float64)

            # update plasmas

    def calculate_updated_radiationfield(self, nubar_estimator, j_estimator):
        """
        Calculate an updated radiation field from the :math:`\\bar{nu}_\\textrm{estimator}` and :math:`\\J_\\textrm{estimator}`
        calculated in the montecarlo simulation. The details of the calculation can be found in the documentation.

        Parameters
        ----------

        nubar_estimator : ~np.ndarray (float)

        j_estimator : ~np.ndarray (float)

        Returns
        -------

        updated_t_rads : ~np.ndarray (float)

        updated_ws : ~np.ndarray (float)

        """
        trad_estimator_constant = 0.260944706 * h / kb # (pi**4 / 15)/ (24*zeta(5))

        updated_t_rads = trad_estimator_constant * nubar_estimator / j_estimator
        updated_ws = j_estimator / (
            4 * constants.cgs.sigma_sb.value * updated_t_rads ** 4 * self.time_of_simulation * self.volumes)

        return updated_t_rads, updated_ws


    def update_plasmas(self, updated_t_rads, updated_ws=None):
        if self.plasma_type == 'lte':
            self.t_rads = updated_t_rads
            for i, (current_plasma, new_trad) in enumerate(zip(self.plasmas, updated_t_rads, updated_ws)):
                current_plasma.set_j_blues()
                current_plasma.update_radiationfield(new_trad)
                self.tau_sobolevs[i] = current_plasma.tau_sobolevs

        elif self.plasma_type == 'nebular':
            self.t_rads = updated_t_rads
            self.ws = updated_ws
            for i, (current_plasma, new_trad, new_ws) in enumerate(zip(self.plasmas, updated_t_rads, updated_ws)):
                current_plasma.set_j_blues()
                current_plasma.update_radiationfield(new_trad, new_ws)
                self.tau_sobolevs[i] = current_plasma.tau_sobolevs


def calculate_atom_number_densities(atom_data, abundances, density):
    """
    Calculates the atom number density, using the following formula, where Z is the atomic number
    and X is the abundance fraction

    .. math::
        N_{Z} = \\frac{\\rho_\\textrm{total}\\times \\textrm{X}_\\textrm{Z}}{m_\\textrm{Z}}

    """

    #Converting abundances



    abundance_fractions = pd.Series(abundances.values(),
        index=pd.Index([atom_data.symbol2atomic_number[item] for item in abundances.keys()],
            dtype=np.int, name='atomic_number'), name='abundance_fraction')



    #Normalizing Abundances

    abundance_sum = abundance_fractions.sum()

    if abs(abundance_sum - 1) > 1e-5:
        logger.warn('Abundances do not add up to 1 (Sum = %.4f). Renormalizing', (abundance_sum))

    abundance_fractions /= abundance_sum

    number_densities = (abundance_fractions * density) /\
                       atom_data.atom_data.ix[abundance_fractions.index]['mass']

    return pd.DataFrame({'abundance_fraction': abundance_fractions, 'number_density': number_densities})


class ModelHistory(object):
    """
    Records the history of the model
    """

    def __init__(self):
        self.plasmas = []

    def store_all(self, radial1d_mdl):
        self.plasmas.append(deepcopy(radial1d_mdl))