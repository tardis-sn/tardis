# building of radial_oned_model

import numpy as np
import plasma, atomic, packet_source
import logging

import pandas as pd
from astropy import constants, units
from copy import deepcopy
import os
import yaml
import pdb

logger = logging.getLogger(__name__)

c = constants.c.cgs.value
h = constants.h.cgs.value
kb = constants.k_B.cgs.value

w_estimator_constant = (c ** 2 / (2 * h)) * (15 / np.pi ** 4) * (h / kb) ** 4 / (4 * np.pi)

synpp_default_yaml = os.path.join(os.path.dirname(__file__), 'data', 'synpp_default.yaml')


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

        self.sigma_thomson = tardis_config.sigma_thomson

        self.spec_nu_bins = np.linspace(tardis_config.spectrum_start_nu, tardis_config.spectrum_end_nu,
                                        tardis_config.spectrum_bins + 1)
        self.spec_nu = self.spec_nu_bins[:-1]

        self.spec_virtual_flux_nu = np.zeros_like(self.spec_nu)


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
            self.ws = np.array([(0.5 * (1 - np.sqrt(1 - self.r_inner[0] ** 2 / self.r_middle[i] ** 2))) if w < 0 \
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
        self.luminosity_inner = 4 * np.pi * constants.sigma_sb.cgs.value * self.r_inner[0] ** 2 * self._t_inner ** 4
        self.time_of_simulation = 1 / self.luminosity_inner


    def create_packets(self, no_of_packets=None):
        #Energy emitted from the inner boundary
        self.emitted_inner_energy = 4 * np.pi * constants.sigma_sb.cgs.value * self.r_inner[0] ** 2 * (
            self.t_inner) ** 4

        if no_of_packets is None:
            no_of_packets = self.no_of_packets
        self.packet_src.create_packets(no_of_packets, self.t_inner)

    def initialize_plasmas(self):
        self.plasmas = []
        self.tau_sobolevs = np.empty((self.no_of_shells, len(self.atom_data.lines)))
        self.line_list_nu = self.atom_data.lines['nu']

        if self.line_interaction_id in (1, 2):
            self.transition_probabilities = []

        if self.plasma_type == 'lte':
            for i, (current_abundances, current_t_rad) in \
                enumerate(zip(self.number_densities, self.t_rads)):
                current_plasma = self.plasma_class(current_abundances, self.atom_data, self.time_explosion,
                                                   nlte_species=self.tardis_config.nlte_species)
                logger.debug('Initializing Shell %d Plasma with T=%.3f' % (i, current_t_rad))
                current_plasma.update_radiationfield(current_t_rad)
                self.tau_sobolevs[i] = current_plasma.tau_sobolevs

                #TODO change this
                if self.line_interaction_id in (1, 2):
                    self.transition_probabilities.append(current_plasma.update_macro_atom().values)

                self.plasmas.append(current_plasma)


        elif self.plasma_type == 'nebular':
            for i, (current_abundances, current_t_rad, current_w) in \
                enumerate(zip(self.number_densities, self.t_rads, self.ws)):
                current_plasma = self.plasma_class(current_abundances, self.atom_data, self.time_explosion,
                                                   nlte_species=self.tardis_config.nlte_species)
                logger.debug('Initializing Shell %d Plasma with T=%.3f W=%.4f' % (i, current_t_rad, current_w))
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
            4 * constants.sigma_sb.cgs.value * updated_t_rads ** 4 * self.time_of_simulation * self.volumes)

        return updated_t_rads, updated_ws


    def update_plasmas(self, updated_t_rads, updated_ws=None):
        if self.plasma_type == 'lte':
            self.t_rads = updated_t_rads
            for i, (current_plasma, new_trad) in enumerate(zip(self.plasmas, updated_t_rads)):
                current_plasma.set_j_blues()

                logger.debug('Updating Shell %d Plasma with T=%.3f' % (i, new_trad))
                current_plasma.update_radiationfield(new_trad)
                self.tau_sobolevs[i] = current_plasma.tau_sobolevs

        elif self.plasma_type == 'nebular':
            self.t_rads = updated_t_rads
            self.ws = updated_ws
            for i, (current_plasma, new_trad, new_ws) in enumerate(zip(self.plasmas, updated_t_rads, updated_ws)):
                current_plasma.set_j_blues()
                logger.debug('Updating Shell %d Plasma with T=%.3f W=%.4f' % (i, new_trad, new_ws))
                current_plasma.update_radiationfield(new_trad, new_ws)
                self.tau_sobolevs[i] = current_plasma.tau_sobolevs

    def calculate_spectrum(self, out_nu, out_energy, distance=None):
        self.out_nu = out_nu
        self.out_energy = out_energy

        if distance is None:
            logger.info('Distance to supernova not selected assuming 10 pc for calculation of spectra')
            distance = units.Quantity(10, 'pc').to('cm').value

        self.spec_flux_nu = np.histogram(out_nu[out_nu > 0], weights=out_energy[out_nu > 0], bins=self.spec_nu_bins)[0]

        flux_scale = self.time_of_simulation * (self.spec_nu[1] - self.spec_nu[0]) * (4 * np.pi * distance ** 2)

        self.spec_flux_nu /= flux_scale

        self.spec_virtual_flux_nu /= flux_scale

        self.spec_reabsorbed_nu = \
            np.histogram(out_nu[out_nu < 0], weights=out_energy[out_nu < 0], bins=self.spec_nu_bins)[0]
        self.spec_reabsorbed_nu /= flux_scale

        self.spec_angstrom = units.Unit('Hz').to('angstrom', self.spec_nu, units.spectral())

        self.spec_flux_angstrom = (self.spec_flux_nu * self.spec_nu ** 2 / constants.c.cgs / 1e8)
        self.spec_reabsorbed_angstrom = (self.spec_reabsorbed_nu * self.spec_nu ** 2 / constants.c.cgs / 1e8)
        self.spec_virtual_flux_angstrom = (self.spec_virtual_flux_nu * self.spec_nu ** 2 / constants.c.cgs / 1e8)

    def create_synpp_yaml(self, fname, lines_db=None):
        if not self.atom_data.has_synpp_refs:
            raise ValueError(
                'The current atom dataset does not contain the necesarry reference files (please contact the authors)')

        self.atom_data.synpp_refs['ref_log_tau'] = -99.0
        for key, value in self.atom_data.synpp_refs.iterrows():
            try:
                tau_sobolev_idx = self.atom_data.lines_index.ix[value['line_id']]
            except KeyError:
                continue

            self.atom_data.synpp_refs['ref_log_tau'].ix[key] = np.log10(self.plasmas[0].tau_sobolevs[tau_sobolev_idx])

        relevant_synpp_refs = self.atom_data.synpp_refs[self.atom_data.synpp_refs['ref_log_tau'] > -50]

        yaml_reference = yaml.load(file(synpp_default_yaml))

        if lines_db is not None:
            yaml_reference['opacity']['line_dir'] = os.path.join(lines_db, 'lines')
            yaml_reference['opacity']['line_dir'] = os.path.join(lines_db, 'refs.dat')

        yaml_setup = yaml_reference['setups'][0]

        yaml_setup['ions'] = []
        yaml_setup['log_tau'] = []
        yaml_setup['active'] = []
        yaml_setup['temp'] = []
        yaml_setup['v_min'] = []
        yaml_setup['v_max'] = []
        yaml_setup['aux'] = []

        for species, synpp_ref in relevant_synpp_refs.iterrows():
            yaml_setup['ions'].append(100 * species[0] + species[1])
            yaml_setup['log_tau'].append(float(synpp_ref['ref_log_tau']))
            yaml_setup['active'].append(True)
            yaml_setup['temp'].append(yaml_setup['t_phot'])
            yaml_setup['v_min'].append(yaml_reference['opacity']['v_ref'])
            yaml_setup['v_max'].append(yaml_reference['grid']['v_outer_max'])
            yaml_setup['aux'].append(1e200)

        yaml.dump(yaml_reference, file(fname, 'w'))


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

    number_densities = (abundance_fractions * density) / \
                       atom_data.atom_data.ix[abundance_fractions.index]['mass']

    return pd.DataFrame({'abundance_fraction': abundance_fractions, 'number_density': number_densities})


class ModelHistory(object):
    """
    Records the history of the model
    """

    def __init__(self, tardis_config):
        self.t_rads = pd.DataFrame(index=np.arange(tardis_config.no_of_shells))
        self.ws = pd.DataFrame(index=np.arange(tardis_config.no_of_shells))
        self.level_populations = {}


    def store_all(self, radial1d_mdl, iteration):
        self.t_rads['iter%d' % iteration] = radial1d_mdl.t_rads
        self.ws['iter%d' % iteration] = radial1d_mdl.ws

        current_level_populations = pd.DataFrame(index=radial1d_mdl.atom_data.levels.index)
        for i, plasma in enumerate(radial1d_mdl.plasmas):
            current_level_populations[i] = plasma.level_populations

        self.level_populations['iter%d' % iteration] = current_level_populations.copy()

    def finalize(self):
        self.level_populations = pd.Panel.from_dict(self.level_populations)


