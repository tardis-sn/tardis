# building of radial_oned_model

import numpy as np
import plasma, packet_source
import logging

import pandas as pd
from pandas import HDFStore
from astropy import constants, units as u
import montecarlo_multizone
import os
import re

import scipy.special

import itertools

logger = logging.getLogger(__name__)

c = constants.c.cgs.value
h = constants.h.cgs.value
kb = constants.k_B.cgs.value

w_estimator_constant = (c ** 2 / (2 * h)) * (15 / np.pi ** 4) * (h / kb) ** 4 / (4 * np.pi)

t_rad_estimator_constant = (np.pi**4 / (15 * 24 * scipy.special.zeta(5, 1))) * h / kb

synpp_default_yaml_fname = os.path.join(os.path.dirname(__file__), 'data', 'synpp_default.yaml')


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

    @classmethod
    def from_h5(cls, buffer_or_fname):
        raise NotImplementedError("This is currently not implemented")


    def __init__(self, tardis_config):
        #final preparation for configuration object
        self.tardis_config = tardis_config
        self.gui = None
        self.converged = False

        self.atom_data = tardis_config.atom_data
        if tardis_config.plasma.type == 'nebular':
            if not self.atom_data.has_zeta_data:
                raise ValueError("Requiring Recombination coefficients Zeta for 'nebular' plasma_type")

        self.packet_src = packet_source.SimplePacketSource.from_wavelength(tardis_config.spectrum.start,
                                                                           tardis_config.spectrum.end)
        self.current_no_of_packets = tardis_config.montecarlo.no_of_packets

        no_of_shells = tardis_config.structure.no_of_shells

        self.t_inner = tardis_config.plasma.t_inner
        self.t_rads = tardis_config.plasma.t_rads

        self.iterations_max_requested = tardis_config.montecarlo.iterations
        self.iterations_remaining = self.iterations_max_requested - 1
        self.iterations_executed = 0


        if tardis_config.montecarlo.convergence.type == 'specific':
            self.global_convergence_parameters = tardis_config.montecarlo.convergence.global_convergence_parameters.config_dict.copy()

        self.t_rads = tardis_config.plasma.t_rads

        self.tau_sobolevs = np.zeros((no_of_shells, len(self.atom_data.lines)))

        self.j_blue_estimators = np.zeros_like(self.tau_sobolevs)
        self.j_blues = np.zeros_like(self.tau_sobolevs)
        j_blues_norm_factor = constants.c.cgs *  tardis_config.supernova.time_explosion / \
                       (4 * np.pi * self.time_of_simulation * tardis_config.structure.volumes.value)
        self.j_blues_norm_factor = j_blues_norm_factor.value.reshape((no_of_shells, 1)) * j_blues_norm_factor.unit

        self.transition_probabilities = np.zeros((no_of_shells, len(self.atom_data.macro_atom_data.lines_idx)))


        self.plasmas = []

        self.ws = (0.5 * (1 - np.sqrt(1 -
                    (tardis_config.structure.r_inner[0] ** 2 / tardis_config.structure.r_middle ** 2).to(1).value)))
        for i in xrange(no_of_shells):
            if tardis_config.plasma.type == 'lte':

                current_plasma_class = plasma.LTEPlasma

            elif tardis_config.plasma.type == 'nebular':

                current_plasma_class = plasma.NebularPlasma

            current_plasma = current_plasma_class(t_rad=self.t_rads[i].value,
                                                     number_density=tardis_config.number_densities.ix[i],
                                                     atom_data=self.atom_data,
                                                     time_explosion=tardis_config.supernova.time_explosion.to('s').value,
                                                     w=self.ws[i], j_blues=self.j_blues[i],
                                                     nlte_config=tardis_config.plasma.nlte, zone_id=i)
            self.plasmas.append(current_plasma)

        self.spectrum = TARDISSpectrum(tardis_config.spectrum.frequency, tardis_config.supernova.distance)
        self.spectrum_virtual = TARDISSpectrum(tardis_config.spectrum.frequency, tardis_config.supernova.distance)
        self.spectrum_reabsorbed = TARDISSpectrum(tardis_config.spectrum.frequency, tardis_config.supernova.distance)



    @property
    def electron_densities(self):
        return np.array([plasma.electron_density for plasma in self.plasmas])

    @property
    def line_interaction_type(self):
        return self._line_interaction_type

    @line_interaction_type.setter
    def line_interaction_type(self, value):
        if value in ['scatter', 'downbranch', 'macroatom']:
            self._line_interaction_type = value
            self.tardis_config.plasma.line_interaction_type = value
            #final preparation for atom_data object - currently building data
            self.atom_data.prepare_atom_data(self.tardis_config.number_densities.columns,
                                             line_interaction_type=self.line_interaction_type, max_ion_number=None,
                                             nlte_species=self.tardis_config.plasma.nlte.species)
        else:
            raise ValueError('line_interaction_type can only be "scatter", "downbranch", or "macroatom"')



    @property
    def t_inner(self):
        return self._t_inner

    @t_inner.setter
    def t_inner(self, value):
        self._t_inner = value
        self.luminosity_inner = (4 * np.pi * constants.sigma_sb.cgs * self.tardis_config.structure.r_inner[0] ** 2 * \
                                self._t_inner ** 4).to('erg/s')
        self.time_of_simulation = (1.0 * u.erg / self.luminosity_inner)


    def calculate_j_blues(self, init_detailed_j_blues=False):
        nus = self.atom_data.lines.nu.values
        t_rads = self.t_rads.value.reshape((self.tardis_config.structure.no_of_shells, 1))
        ws = self.ws.reshape((self.tardis_config.structure.no_of_shells, 1))
        radiative_rates_type = self.tardis_config.plasma.radiative_rates_type
        w_epsilon = self.tardis_config.plasma.w_epsilon

        if radiative_rates_type == 'lte':
            logger.info('Calculating J_blues for radiate_rates_type=lte')
            self.j_blues = plasma.intensity_black_body(nus, t_rads)
        elif radiative_rates_type == 'nebular' or init_detailed_j_blues:
            logger.info('Calculating J_blues for radiate_rates_type=nebular')
            self.j_blues = ws * plasma.intensity_black_body(nus, t_rads)
        elif radiative_rates_type == 'detailed':
            logger.info('Calculating J_blues for radiate_rates_type=nebular')
            self.j_blues = self.j_blue_estimators * self.j_blues_norm_factor.value
            for i in xrange(self.tardis_config.structure.no_of_shells):
                zero_j_blues = self.j_blues[i] == 0.0
                self.j_blues[i][zero_j_blues] = w_epsilon * plasma.intensity_black_body(
                    self.atom_data.lines.nu.values[zero_j_blues], self.t_rads.value[i])





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


        updated_t_rads = t_rad_estimator_constant * nubar_estimator / j_estimator
        updated_ws = j_estimator / (
            4 * constants.sigma_sb.cgs.value * updated_t_rads ** 4 * self.time_of_simulation.value
            * self.tardis_config.structure.volumes.value)

        return updated_t_rads * u.K, updated_ws



    def update_plasmas(self):
        for i in xrange(self.tardis_config.structure.no_of_shells):
        #for i, (current_plasma, new_trad, new_ws) in enumerate(zip(self.plasmas, self.t_rads, self.ws)):
            logger.debug('Updating Shell %d Plasma with T=%.3f W=%.4f' % (i, self.t_rads[i].value, self.ws[i]))
            self.plasmas[i].set_j_blues(self.j_blues[i])
            if self.tardis_config.plasma.type == 'lte':
                current_ws = 1.0
            elif self.tardis_config.plasma.type == 'nebular':
                current_ws = self.ws[i]

            self.plasmas[i].update_radiationfield(self.t_rads[i].value, w=current_ws)
            self.tau_sobolevs[i] = self.plasmas[i].tau_sobolevs

            if self.tardis_config.plasma.line_interaction_type in ('downbranch', 'macroatom'):
                self.transition_probabilities[i] = self.plasmas[i].calculate_transition_probabilities()




    def simulate(self, update_radiation_field=True, enable_virtual=False, initialize_j_blues=False):
        """
        Run a simulation
        """

        self.packet_src.create_packets(self.current_no_of_packets, self.t_inner.value)
        self.calculate_j_blues(init_detailed_j_blues=initialize_j_blues)
        self.update_plasmas()


        if update_radiation_field:
            self.update_radiationfield()


        if enable_virtual:
            no_of_virtual_packets = self.tardis_config.montecarlo.no_of_virtual_packets
        else:
            no_of_virtual_packets = 0
        if np.any(np.isnan(self.tau_sobolevs)) or np.any(np.isinf(self.tau_sobolevs)) or np.any(np.isneginf(self.tau_sobolevs)):
            raise ValueError('Some values are nan, inf, -inf in tau_sobolevs. Something went wrong!')


        self.montecarlo_virtual_luminosity = np.zeros_like(self.spectrum.frequency.value)
        montecarlo_nu, montecarlo_energies, self.j_estimators, self.nubar_estimators, \
        last_line_interaction_in_id, last_line_interaction_out_id, \
        self.last_interaction_type, self.last_line_interaction_shell_id = \
            montecarlo_multizone.montecarlo_radial1d(self,
                                                     virtual_packet_flag=no_of_virtual_packets)

        self.montecarlo_nu = montecarlo_nu * u.Hz
        self.montecarlo_luminosity = montecarlo_energies *  1 * u.erg / self.time_of_simulation

        montecarlo_reabsorbed_luminosity = -np.histogram(self.montecarlo_nu.value[self.montecarlo_luminosity.value < 0],
                                         weights=self.montecarlo_luminosity.value[self.montecarlo_luminosity.value < 0],
                                         bins=self.tardis_config.spectrum.frequency.value)[0] \
                                      * self.montecarlo_luminosity.unit

        montecarlo_emitted_luminosity = np.histogram(self.montecarlo_nu.value[self.montecarlo_luminosity.value >= 0],
                                         weights=self.montecarlo_luminosity.value[self.montecarlo_luminosity.value >= 0],
                                         bins=self.tardis_config.spectrum.frequency.value)[0] \
                                   * self.montecarlo_luminosity.unit



        self.spectrum.update_luminosity(montecarlo_emitted_luminosity)
        self.spectrum_reabsorbed.update_luminosity(montecarlo_reabsorbed_luminosity)


        if no_of_virtual_packets > 0:
            self.montecarlo_virtual_luminosity = self.montecarlo_virtual_luminosity \
                                                 * 1 * u.erg / self.time_of_simulation
            self.spectrum_virtual.update_luminosity(self.montecarlo_virtual_luminosity)



        self.last_line_interaction_in_id = self.atom_data.lines_index.index.values[last_line_interaction_in_id]
        self.last_line_interaction_in_id[last_line_interaction_in_id == -1] = -1
        self.last_line_interaction_out_id = self.atom_data.lines_index.index.values[last_line_interaction_out_id]
        self.last_line_interaction_out_id[last_line_interaction_out_id == -1] = -1

        self.iterations_executed += 1
        self.iterations_remaining -= 1

        if self.gui is not None:
            self.gui.update_data(self)
            self.gui.show()



    def update_radiationfield(self, log_sampling=5):
        """
        Updating radiation field
        """

        updated_t_rads, updated_ws = self.calculate_updated_radiationfield(self.nubar_estimators, self.j_estimators)
        old_t_rads = self.t_rads.copy()
        old_ws = self.ws.copy()
        old_t_inner = self.t_inner
        luminosity_wavelength_filter = (self.montecarlo_nu > self.tardis_config.supernova.luminosity_nu_start) & \
                            (self.montecarlo_nu < self.tardis_config.supernova.luminosity_nu_end)
        emitted_luminosity = np.sum(self.montecarlo_luminosity.value[(self.montecarlo_luminosity.value >= 0) &
                                                         luminosity_wavelength_filter]) \
                             * self.montecarlo_luminosity.unit

        absorbed_luminosity = -np.sum(self.montecarlo_luminosity.value[(self.montecarlo_luminosity.value < 0) &
                                                          luminosity_wavelength_filter]) \
                              * self.montecarlo_luminosity.unit
        updated_t_inner = self.t_inner \
                          * (emitted_luminosity / self.tardis_config.supernova.luminosity_requested).to(1).value ** -.25

        convergence_t_rads = (abs(old_t_rads - updated_t_rads) / updated_t_rads).value
        convergence_ws = (abs(old_ws - updated_ws) / updated_ws)
        convergence_t_inner = (abs(old_t_inner - updated_t_inner) / updated_t_inner).value

        convergence_section = self.tardis_config.montecarlo.convergence
        if convergence_section.type == 'damped' or convergence_section.type == 'specific':
            self.t_rads += convergence_section.t_rad.damping_constant * (updated_t_rads - self.t_rads)
            self.ws += convergence_section.w.damping_constant * (updated_ws - self.ws)
            self.t_inner += convergence_section.w.damping_constant * (updated_t_inner - self.t_inner)

        if convergence_section.type == 'specific':

            t_rad_converged = (float(np.sum(convergence_t_rads < convergence_section.t_rad['threshold'])) \
                               / self.tardis_config.structure.no_of_shells) >= convergence_section.t_rad['fraction']

            w_converged = (float(np.sum(convergence_t_rads < convergence_section.w['threshold'])) \
                           / self.tardis_config.structure.no_of_shells) >= convergence_section.w['fraction']

            t_inner_converged = convergence_t_inner < convergence_section.t_inner['threshold']

            if t_rad_converged and t_inner_converged and w_converged:
                if not self.converged:
                    self.converged = True
                    self.iterations_remaining = self.global_convergence_parameters['hold']

            else:
                if self.converged:
                    self.iterations_remaining = self.iterations_max_requested - self.iterations_executed
                    self.converged = False

        self.temperature_logging = pd.DataFrame(
            {'t_rads': old_t_rads.value, 'updated_t_rads': updated_t_rads.value,
             'converged_t_rads': convergence_t_rads, 'new_trads': self.t_rads.value, 'ws': old_ws,
             'updated_ws': updated_ws, 'converged_ws': convergence_ws,
             'new_ws': self.ws})

        self.temperature_logging.index.name = 'Shell'

        temperature_logging = str(self.temperature_logging[::log_sampling])

        temperature_logging = ''.join(['\t%s\n' % item for item in temperature_logging.split('\n')])

        logger.info('Plasma stratification:\n%s\n', temperature_logging)
        logger.info("Luminosity emitted = %.5e Luminosity absorbed = %.5e Luminosity requested = %.5e",
                    emitted_luminosity.value, absorbed_luminosity.value,
                    self.tardis_config.supernova.luminosity_requested.value)
        logger.info('Calculating new t_inner = %.3f', updated_t_inner.value)

    def save_spectra(self, fname):
        self.spectrum.to_ascii(fname)
        self.spectrum_virtual.to_ascii('virtual_' + fname)


    def to_hdf5(self, buffer_or_fname, path='', close_h5=True):
        if isinstance(buffer_or_fname, basestring):
            hdf_store = pd.HDFStore(buffer_or_fname)
        elif isinstance(buffer_or_fname, pd.HDFStore):
            hdf_store = buffer_or_fname
        else:
            raise IOError('Please specify either a filename or an HDFStore')
        logger.info('Writing to path %s', path)


        level_populations_path = os.path.join(path, 'level_populations')
        level_populations = np.zeros((self.tardis_config.structure.no_of_shells,
                                      len(self.plasmas[0].level_populations)))

        ion_populations_path = os.path.join(path, 'ion_populations')
        ion_populations = np.zeros((self.tardis_config.structure.no_of_shells,
                                    len(self.plasmas[0].ion_populations)))

        for i in xrange(self.tardis_config.structure.no_of_shells):
            level_populations[i] = self.plasmas[i].level_populations.values
            ion_populations[i] = self.plasmas[i].ion_populations.values



        pd.DataFrame(level_populations.transpose(), index=self.atom_data.levels.index).to_hdf(hdf_store,
                                                                                              level_populations_path)
        pd.DataFrame(ion_populations.transpose(), index=self.plasmas[0].ion_populations.index).to_hdf(hdf_store,
                                                                                                      ion_populations_path)

        tau_sobolevs_path = os.path.join(path, 'tau_sobolevs')
        pd.DataFrame(self.tau_sobolevs.transpose(), index=self.atom_data.lines.index).to_hdf(hdf_store,
                                                                                             tau_sobolevs_path)

        j_blues_path = os.path.join(path, 'j_blues')
        pd.DataFrame(self.j_blues.transpose(), index=self.atom_data.lines.index).to_hdf(hdf_store, j_blues_path)

        t_rads_path = os.path.join(path, 't_rads')
        pd.Series(self.t_rads.value).to_hdf(hdf_store, t_rads_path)

        ws_path = os.path.join(path, 'ws')
        pd.Series(self.ws).to_hdf(hdf_store, ws_path)

        electron_densities_path = os.path.join(path, 'electron_densities')
        pd.Series(self.electron_densities).to_hdf(hdf_store, electron_densities_path)

        last_line_interaction_in_id_path = os.path.join(path, 'last_line_interaction_in_id')
        pd.Series(self.last_line_interaction_in_id).to_hdf(hdf_store, last_line_interaction_in_id_path)

        last_line_interaction_out_id_path = os.path.join(path, 'last_line_interaction_out_id')
        pd.Series(self.last_line_interaction_out_id).to_hdf(hdf_store, last_line_interaction_out_id_path)

        last_line_interaction_shell_id_path = os.path.join(path, 'last_line_interaction_shell_id')
        pd.Series(self.last_line_interaction_shell_id).to_hdf(hdf_store, last_line_interaction_shell_id_path)

        luminosity_density = pd.DataFrame.from_dict(dict(wave=self.spectrum.wavelength.value,
                                                         flux=self.spectrum.luminosity_density_lambda.value))

        luminosity_density.to_hdf(hdf_store, os.path.join(path, 'luminosity_density'))

        if self.spectrum_virtual.luminosity_density_lambda is not None:
            luminosity_density_virtual = pd.DataFrame.from_dict(dict(wave=self.spectrum_virtual.wavelength.value,
                                                           flux=self.spectrum_virtual.luminosity_density_lambda.value))
            luminosity_density_virtual.to_hdf(hdf_store, os.path.join(path, 'luminosity_density_virtual'))

        hdf_store.flush()
        if close_h5:
            hdf_store.close()
        else:
            return hdf_store




class TARDISHistory(object):
    """
    Records the history of the model
    """


    @classmethod
    def from_hdf5(cls, fname):
        hdf_store = HDFStore(fname)
        iterations = []
        for key in hdf_store.keys():
            iterations.append(int(re.match('model(\d+)', key.split('/')[1]).groups()[0]))

        iterations = np.sort(np.unique(iterations))

        history = cls()
        no_of_shell = len(hdf_store['model1/t_rads'].index)
        t_rads_dict = {}
        ws_dict = {}
        level_populations_dict = {}
        ion_populations_dict = {}

        history.iterations = iterations

        for iter in iterations:
            current_iter = 'iter%d' % iter
            t_rads_dict[current_iter] = hdf_store['model%d/t_rads' % iter]
            ws_dict[current_iter] = hdf_store['model%d/ws' % iter]
            level_populations_dict[current_iter] = hdf_store['model%d/level_populations' % iter]
            ion_populations_dict[current_iter] = hdf_store['model%d/ion_populations' % iter]

            for index in ion_populations_dict[current_iter].index:
                level_populations_dict[current_iter].ix[index].update(level_populations_dict[current_iter].ix[index] /
                                                                      ion_populations_dict[current_iter].ix[index])





        history.t_rads = pd.DataFrame(t_rads_dict)
        history.ws = pd.DataFrame(ws_dict)
        history.level_populations = pd.Panel(level_populations_dict)
        history.ion_populations = pd.Panel(ion_populations_dict)

        return history

    def plot_level_evolution(self, level_index, shell):
        for iter in self.iterations:
            pass





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
