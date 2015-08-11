# This module contains the model class

import logging
import os
import itertools

import numpy as np
import pandas as pd
from astropy import constants, units as u
import scipy.special

from util import intensity_black_body
from tardis import packet_source, plasma_array
from tardis.montecarlo import montecarlo
from tardis.montecarlo.base import MontecarloRunner



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
        self.gui = None
        self.converged = False
        self.atom_data = tardis_config.atom_data
        selected_atomic_numbers = self.tardis_config.abundances.index
        self.atom_data.prepare_atom_data(selected_atomic_numbers,
                                         line_interaction_type=tardis_config.plasma.line_interaction_type,
                                         nlte_species=tardis_config.plasma.nlte.species)

        if tardis_config.plasma.ionization == 'nebular':
            if not self.atom_data.has_zeta_data:
                raise ValueError("Requiring Recombination coefficients Zeta for 'nebular' plasma ionization")

        self.packet_src = packet_source.SimplePacketSource.from_wavelength(tardis_config.montecarlo.black_body_sampling.start,
                                                                           tardis_config.montecarlo.black_body_sampling.end,
                                                                           blackbody_sampling=tardis_config.montecarlo.black_body_sampling.samples,
                                                                           seed=self.tardis_config.montecarlo.seed)
        self.current_no_of_packets = tardis_config.montecarlo.no_of_packets

        self.t_inner = tardis_config.plasma.t_inner
        self.t_rads = tardis_config.plasma.t_rads

        self.iterations_max_requested = tardis_config.montecarlo.iterations
        self.iterations_remaining = self.iterations_max_requested
        self.iterations_executed = 0


        if tardis_config.montecarlo.convergence_strategy.type == 'specific':
            self.global_convergence_parameters = (tardis_config.montecarlo.
                                                  convergence_strategy.
                                                  deepcopy())

        self.t_rads = tardis_config.plasma.t_rads
        t_inner_lock_cycle = [False] * (tardis_config.montecarlo.
                                        convergence_strategy.
                                        lock_t_inner_cycles)
        t_inner_lock_cycle[0] = True
        self.t_inner_update = itertools.cycle(t_inner_lock_cycle)



        self.ws = (0.5 * (1 - np.sqrt(1 -
                    (tardis_config.structure.r_inner[0] ** 2 / tardis_config.structure.r_middle ** 2).to(1).value)))


        self.plasma_array = plasma_array.BasePlasmaArray(tardis_config.number_densities, tardis_config.atom_data,
                                                         tardis_config.supernova.time_explosion.to('s').value,
                                                         nlte_config=tardis_config.plasma.nlte,
                                                         delta_treatment=tardis_config.plasma.delta_treatment,
                                                         ionization_mode=tardis_config.plasma.ionization,
                                                         excitation_mode=tardis_config.plasma.excitation)




        self.spectrum = TARDISSpectrum(tardis_config.spectrum.frequency, tardis_config.supernova.distance)
        self.spectrum_virtual = TARDISSpectrum(tardis_config.spectrum.frequency, tardis_config.supernova.distance)
        self.spectrum_reabsorbed = TARDISSpectrum(tardis_config.spectrum.frequency, tardis_config.supernova.distance)
        self.runner = MontecarloRunner()




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
                                self.t_inner ** 4).to('erg/s')
        self.time_of_simulation = (1.0 * u.erg / self.luminosity_inner)
        self.j_blues_norm_factor = constants.c.cgs *  self.tardis_config.supernova.time_explosion / \
                       (4 * np.pi * self.time_of_simulation * self.tardis_config.structure.volumes)


    def calculate_j_blues(self, init_detailed_j_blues=False):
        nus = self.atom_data.lines.nu.values
        radiative_rates_type = self.tardis_config.plasma.radiative_rates_type
        w_epsilon = self.tardis_config.plasma.w_epsilon

        if radiative_rates_type == 'lte':
            logger.info('Calculating J_blues for radiative_rates_type=lte')
            j_blues = intensity_black_body(nus[np.newaxis].T, self.t_rads.value)
            self.j_blues = pd.DataFrame(j_blues, index=self.atom_data.lines.index, columns=np.arange(len(self.t_rads)))
        elif radiative_rates_type == 'dilute-blackbody' or init_detailed_j_blues:
            logger.info('Calculating J_blues for radiative_rates_type=dilute-blackbody')
            j_blues = self.ws * intensity_black_body(nus[np.newaxis].T, self.t_rads.value)
            self.j_blues = pd.DataFrame(j_blues, index=self.atom_data.lines.index, columns=np.arange(len(self.t_rads)))

        elif radiative_rates_type == 'detailed':
            logger.info('Calculating J_blues for radiate_rates_type=detailed')

            self.j_blues = pd.DataFrame(self.j_blue_estimators.transpose() * self.j_blues_norm_factor.value,
                                        index=self.atom_data.lines.index, columns=np.arange(len(self.t_rads)))
            for i in xrange(self.tardis_config.structure.no_of_shells):
                zero_j_blues = self.j_blues[i] == 0.0
                self.j_blues[i][zero_j_blues] = w_epsilon * intensity_black_body(
                    self.atom_data.lines.nu.values[zero_j_blues], self.t_rads.value[i])

        else:
            raise ValueError('radiative_rates_type type unknown - %s', radiative_rates_type)

    def update_plasmas(self, initialize_nlte=False):

        self.plasma_array.update_radiationfield(self.t_rads.value, self.ws, j_blues=self.j_blues,
                                        initialize_nlte=initialize_nlte)


        if self.tardis_config.plasma.line_interaction_type in ('downbranch', 'macroatom'):
            self.transition_probabilities = self.plasma_array.calculate_transition_probabilities()


    def update_radiationfield(self, log_sampling=5):
        """
        Updating radiation field
        """
        convergence_section = self.tardis_config.montecarlo.convergence_strategy
        updated_t_rads, updated_ws = (
            self.runner.calculate_radiationfield_properties())
        old_t_rads = self.t_rads.copy()
        old_ws = self.ws.copy()
        old_t_inner = self.t_inner
        luminosity_wavelength_filter = (self.montecarlo_nu > self.tardis_config.supernova.luminosity_nu_start) & \
                            (self.montecarlo_nu < self.tardis_config.supernova.luminosity_nu_end)
        emitted_filter = self.montecarlo_luminosity.value >= 0
        emitted_luminosity = np.sum(self.montecarlo_luminosity.value[emitted_filter & luminosity_wavelength_filter]) \
                             * self.montecarlo_luminosity.unit

        absorbed_luminosity = -np.sum(self.montecarlo_luminosity.value[~emitted_filter & luminosity_wavelength_filter]) \
                              * self.montecarlo_luminosity.unit
        updated_t_inner = self.t_inner \
                          * (emitted_luminosity / self.tardis_config.supernova.luminosity_requested).to(1).value \
                            ** convergence_section.t_inner_update_exponent

        #updated_t_inner = np.max([np.min([updated_t_inner, 30000]), 3000])

        convergence_t_rads = (abs(old_t_rads - updated_t_rads) / updated_t_rads).value
        convergence_ws = (abs(old_ws - updated_ws) / updated_ws)
        convergence_t_inner = (abs(old_t_inner - updated_t_inner) / updated_t_inner).value


        if convergence_section.type == 'damped' or convergence_section.type == 'specific':
            self.t_rads += convergence_section.t_rad.damping_constant * (updated_t_rads - self.t_rads)
            self.ws += convergence_section.w.damping_constant * (updated_ws - self.ws)
            if self.t_inner_update.next():
                t_inner_new = self.t_inner + convergence_section.t_inner.damping_constant * (updated_t_inner - self.t_inner)
            else:
                t_inner_new = self.t_inner


        if convergence_section.type == 'specific':

            t_rad_converged = (float(np.sum(convergence_t_rads < convergence_section.t_rad['threshold'])) \
                               / self.tardis_config.structure.no_of_shells) >= convergence_section.t_rad['fraction']

            w_converged = (float(np.sum(convergence_t_rads < convergence_section.w['threshold'])) \
                           / self.tardis_config.structure.no_of_shells) >= convergence_section.w['fraction']

            t_inner_converged = convergence_t_inner < convergence_section.t_inner['threshold']

            if t_rad_converged and t_inner_converged and w_converged:
                if not self.converged:
                    self.converged = True
                    self.iterations_remaining = self.global_convergence_parameters['hold_iterations']

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


        return t_inner_new


    def simulate(self, update_radiation_field=True, enable_virtual=False, initialize_j_blues=False,
                 initialize_nlte=False):
        """
        Run a simulation
        """

        if update_radiation_field:
            t_inner_new = self.update_radiationfield()
        else:
            t_inner_new = self.t_inner

        self.calculate_j_blues(init_detailed_j_blues=initialize_j_blues)
        self.update_plasmas(initialize_nlte=initialize_nlte)


        self.t_inner = t_inner_new

        self.packet_src.create_packets(self.current_no_of_packets, self.t_inner.value)

        if enable_virtual:
            no_of_virtual_packets = self.tardis_config.montecarlo.no_of_virtual_packets
        else:
            no_of_virtual_packets = 0
        if np.any(np.isnan(self.plasma_array.tau_sobolevs.values)) or np.any(np.isinf(self.plasma_array.tau_sobolevs.values)) \
            or np.any(np.isneginf(self.plasma_array.tau_sobolevs.values)):
            raise ValueError('Some tau_sobolevs are nan, inf, -inf in tau_sobolevs. Something went wrong!')

        self.j_blue_estimators = np.zeros((len(self.t_rads), len(self.atom_data.lines)))
        self.montecarlo_virtual_luminosity = np.zeros_like(self.spectrum.frequency.value)

        self.runner.run(self, no_of_virtual_packets=no_of_virtual_packets,
                        nthreads=self.tardis_config.montecarlo.nthreads) #self = model


        (montecarlo_nu, montecarlo_energies, self.j_estimators,
         self.nubar_estimators, last_line_interaction_in_id,
         last_line_interaction_out_id, self.last_interaction_type,
         self.last_line_interaction_shell_id) = self.runner.legacy_return()

        if np.sum(montecarlo_energies < 0) == len(montecarlo_energies):
            logger.critical("No r-packet escaped through the outer boundary.")

        self.montecarlo_nu = self.runner.packet_nu
        self.montecarlo_luminosity = self.runner.packet_luminosity



        montecarlo_reabsorbed_luminosity = np.histogram(
            self.runner.reabsorbed_packet_nu,
            weights=self.runner.reabsorbed_packet_luminosity,
            bins=self.tardis_config.spectrum.frequency.value)[0] * u.erg / u.s



        montecarlo_emitted_luminosity = np.histogram(
            self.runner.emitted_packet_nu,
            weights=self.runner.emitted_packet_luminosity,
            bins=self.tardis_config.spectrum.frequency.value)[0] * u.erg / u.s



        self.spectrum.update_luminosity(montecarlo_emitted_luminosity)
        self.spectrum_reabsorbed.update_luminosity(montecarlo_reabsorbed_luminosity)


        if no_of_virtual_packets > 0:
            self.montecarlo_virtual_luminosity = self.montecarlo_virtual_luminosity \
                                                 * 1 * u.erg / self.time_of_simulation
            self.spectrum_virtual.update_luminosity(self.montecarlo_virtual_luminosity)



        self.last_line_interaction_in_id = self.atom_data.lines_index.index.values[last_line_interaction_in_id]
        self.last_line_interaction_in_id = self.last_line_interaction_in_id[last_line_interaction_in_id != -1]
        self.last_line_interaction_out_id = self.atom_data.lines_index.index.values[last_line_interaction_out_id]
        self.last_line_interaction_out_id = self.last_line_interaction_out_id[last_line_interaction_out_id != -1]
        self.last_line_interaction_angstrom = self.montecarlo_nu[last_line_interaction_in_id != -1].to('angstrom',
                                                                                                       u.spectral())


        self.iterations_executed += 1
        self.iterations_remaining -= 1

        if self.gui is not None:
            self.gui.update_data(self)
            self.gui.show()



    def save_spectra(self, fname):
        self.spectrum.to_ascii(fname)
        self.spectrum_virtual.to_ascii('virtual_' + fname)


    def to_hdf5(self, buffer_or_fname, path='', close_h5=True):
        """
            This allows the model to be written to an HDF5 file for later analysis. Currently, the saved properties
            are specified hard coded in include_from_model_in_hdf5. This is a dict where the key corresponds to the
            name of the property and the value describes the type. If the value is None the property can be dumped
            to hdf via its attribute to_hdf or by converting it to a pd.DataFrame. For more complex properties
            which can not simply be dumped to an hdf file the dict can contain a function which is called with
            the parameters key, path, and  hdf_store. This function then should dump the data to the given
            hdf_store object. To dump  properties of sub-properties of  the model, you can use a dict as value.
            This dict is then treated in the same way as described above.

        Parameters
        ----------

        buffer_or_fname: buffer or ~str
            buffer or filename for HDF5 file (see pandas.HDFStore for description)
        path: ~str, optional
            path in the HDF5 file
        close_h5: ~bool
            close the HDF5 file or not.
        """


        # Functions to save properties of the model without to_hdf attribute and no simple conversion to a pd.DataFrame.
        #This functions are always called with the parameters key, path and,  hdf_store.
        def _save_luminosity_density(key, path, hdf_store):

            luminosity_density = pd.DataFrame.from_dict(dict(wave=self.spectrum.wavelength.value,
                                                             flux=self.spectrum.luminosity_density_lambda.value))
            luminosity_density.to_hdf(hdf_store, os.path.join(path, key))

        def _save_spectrum_virtual(key, path, hdf_store):
            if self.spectrum_virtual.luminosity_density_lambda is not None:
                luminosity_density_virtual = pd.DataFrame.from_dict(dict(wave=self.spectrum_virtual.wavelength.value,
                                                                         flux=self.spectrum_virtual.luminosity_density_lambda.value))
                luminosity_density_virtual.to_hdf(hdf_store, os.path.join(path, key))

        def _save_configuration_dict(key, path, hdf_store):
            configuration_dict = dict(t_inner=self.t_inner.value)
            configuration_dict_path = os.path.join(path, 'configuration')
            pd.Series(configuration_dict).to_hdf(hdf_store, configuration_dict_path)

        include_from_plasma_ = {'level_populations': None, 'ion_populations': None, 'tau_sobolevs': None,
                                'electron_densities': None,
                                't_rads': None, 'ws': None}
        include_from_model_in_hdf5 = {'plasma_array': include_from_plasma_, 'j_blues': None,
                                      'last_line_interaction_in_id': None,
                                      'last_line_interaction_out_id': None,
                                      'last_line_interaction_shell_id': None, 'montecarlo_nu': None,
                                      'luminosity_density': _save_luminosity_density,
                                      'luminosity_density_virtual': _save_spectrum_virtual,
                                      'configuration_dict': _save_configuration_dict,
                                      'last_line_interaction_angstrom': None}

        if isinstance(buffer_or_fname, basestring):
            hdf_store = pd.HDFStore(buffer_or_fname)
        elif isinstance(buffer_or_fname, pd.HDFStore):
            hdf_store = buffer_or_fname
        else:
            raise IOError('Please specify either a filename or an HDFStore')
        logger.info('Writing to path %s', path)

        def _get_hdf5_path(path, property_name):
            return os.path.join(path, property_name)

        def _to_smallest_pandas(object):
            try:
                return pd.Series(object)
            except Exception:
                return pd.DataFrame(object)


        def _save_model_property(object, property_name, path, hdf_store):
            property_path = _get_hdf5_path(path, property_name)

            try:
                object.to_hdf(hdf_store, property_path)
            except AttributeError:
                _to_smallest_pandas(object).to_hdf(hdf_store, property_path)


        for key in include_from_model_in_hdf5:
            if include_from_model_in_hdf5[key] is None:
                _save_model_property(getattr(self, key), key, path, hdf_store)
            elif callable(include_from_model_in_hdf5[key]):
                include_from_model_in_hdf5[key](key, path, hdf_store)
            else:
                try:
                    for subkey in include_from_model_in_hdf5[key]:
                        if include_from_model_in_hdf5[key][subkey] is None:
                            _save_model_property(getattr(getattr(self, key), subkey), subkey, os.path.join(path, key),
                                                 hdf_store)
                        elif callable(include_from_model_in_hdf5[key][subkey]):
                            include_from_model_in_hdf5[key][subkey](subkey, os.path.join(path, key), hdf_store)
                        else:
                            logger.critical('Can not save %s', str(os.path.join(path, key, subkey)))
                except:
                    logger.critical('An error occurred while dumping %s to HDF.', str(os.path.join(path, key)))


        hdf_store.flush()
        if close_h5:
            hdf_store.close()
        else:
            return hdf_store

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
