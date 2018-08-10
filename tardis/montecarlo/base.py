import os
import logging
import warnings

from astropy import units as u, constants as const

from scipy.special import zeta
from spectrum import TARDISSpectrum

from tardis.util.base import quantity_linspace
from tardis.io.util import HDFWriterMixin
from tardis.montecarlo import montecarlo, packet_source
from tardis.montecarlo.formal_integral import FormalIntegrator

import numpy as np

logger = logging.getLogger(__name__)


class MontecarloRunner(HDFWriterMixin):
    """
    This class is designed as an interface between the Python part and the
    montecarlo C-part
    """
    hdf_properties = ['output_nu', 'output_energy', 'nu_bar_estimator',
                      'j_estimator', 'montecarlo_virtual_luminosity',
                      'last_interaction_in_nu',
                      'last_interaction_type',
                      'last_line_interaction_in_id',
                      'last_line_interaction_out_id',
                      'last_line_interaction_shell_id',
                      'packet_luminosity', 'spectrum',
                      'spectrum_virtual', 'spectrum_reabsorbed']
    hdf_name = 'runner'
    w_estimator_constant = ((const.c ** 2 / (2 * const.h)) *
                            (15 / np.pi ** 4) * (const.h / const.k_B) ** 4 /
                            (4 * np.pi)).cgs.value

    t_rad_estimator_constant = ((np.pi**4 / (15 * 24 * zeta(5, 1))) *
                                (const.h / const.k_B)).cgs.value

    def __init__(self, seed, spectrum_frequency, virtual_spectrum_range,
                 sigma_thomson, enable_reflective_inner_boundary,
                 inner_boundary_albedo, line_interaction_type):

        self.seed = seed
        self.packet_source = packet_source.BlackBodySimpleSource(seed)
        self.spectrum_frequency = spectrum_frequency
        self.virtual_spectrum_range = virtual_spectrum_range
        self.sigma_thomson = sigma_thomson
        self.enable_reflective_inner_boundary = enable_reflective_inner_boundary
        self.inner_boundary_albedo = inner_boundary_albedo
        self.line_interaction_type = line_interaction_type
        self._integrator = None
        self._spectrum_integrated = None

    def _initialize_estimator_arrays(self, no_of_shells, tau_sobolev_shape):
        """
        Initialize the output arrays of the montecarlo simulation.

        Parameters
        ----------

        model: ~Radial1DModel
        """

        # Estimators
        self.j_estimator = np.zeros(no_of_shells, dtype=np.float64)
        self.nu_bar_estimator = np.zeros(no_of_shells, dtype=np.float64)
        self.j_blue_estimator = np.zeros(tau_sobolev_shape)
        self.Edotlu_estimator = np.zeros(tau_sobolev_shape)

    def _initialize_geometry_arrays(self, model):
        """
        Generate the cgs like geometry arrays for the montecarlo part

        Parameters
        ----------

        model : model.Radial1DModel
        """
        self.r_inner_cgs = model.r_inner.to('cm').value
        self.r_outer_cgs = model.r_outer.to('cm').value
        self.v_inner_cgs = model.v_inner.to('cm/s').value

    def _initialize_packets(self, T, no_of_packets):
        nus, mus, energies = self.packet_source.create_packets(
                T,
                no_of_packets
                )
        self.input_nu = nus
        self.input_mu = mus
        self.input_energy = energies

        self._output_nu = np.ones(no_of_packets, dtype=np.float64) * -99.0
        self._output_energy = np.ones(no_of_packets, dtype=np.float64) * -99.0

        self.last_line_interaction_in_id = -1 * np.ones(
            no_of_packets, dtype=np.int64)
        self.last_line_interaction_out_id = -1 * np.ones(
            no_of_packets, dtype=np.int64)
        self.last_line_interaction_shell_id = -1 * np.ones(
            no_of_packets, dtype=np.int64)
        self.last_interaction_type = -1 * np.ones(
                no_of_packets, dtype=np.int64)
        self.last_interaction_in_nu = np.zeros(no_of_packets, dtype=np.float64)

        self._montecarlo_virtual_luminosity = u.Quantity(
                np.zeros_like(self.spectrum_frequency.value),
                'erg / s'
                )

    @property
    def spectrum(self):
        return TARDISSpectrum(
                self.spectrum_frequency,
                self.montecarlo_emitted_luminosity)

    @property
    def spectrum_reabsorbed(self):
        return TARDISSpectrum(
                self.spectrum_frequency,
                self.montecarlo_reabsorbed_luminosity)

    @property
    def spectrum_virtual(self):
        if np.all(self.montecarlo_virtual_luminosity == 0):
            warnings.warn(
                    "MontecarloRunner.spectrum_virtual"
                    "is zero. Please run the montecarlo simulation with"
                    "no_of_virtual_packets > 0", UserWarning)

        return TARDISSpectrum(
                self.spectrum_frequency,
                self.montecarlo_virtual_luminosity)

    @property
    def spectrum_integrated(self):
        if self._spectrum_integrated is None:
            self._spectrum_integrated = self.integrator.calculate_spectrum(
                self.spectrum_frequency[:-1])
        return self._spectrum_integrated

    @property
    def integrator(self):
        if self._integrator is None:
            warnings.warn(
                    "MontecarloRunner.integrator: "
                    "The FormalIntegrator is not yet available."
                    "Please run the montecarlo simulation at least once.",
                    UserWarning)
        return self._integrator

    def run(self, model, plasma, no_of_packets,
            no_of_virtual_packets=0, nthreads=1,
            last_run=False):
        """
        Run the montecarlo calculation

        Parameters
        ----------
        model : tardis.model.Radial1DModel
        plasma : tardis.plasma.BasePlasma
        no_of_packets : int
        no_of_virtual_packets : int
        nthreads : int
        last_run : bool

        Returns
        -------
        None
        """
        self._integrator = FormalIntegrator(
                model,
                plasma,
                self)
        self.time_of_simulation = self.calculate_time_of_simulation(model)
        self.volume = model.volume
        self._initialize_estimator_arrays(self.volume.shape[0],
                                          plasma.tau_sobolevs.shape)
        self._initialize_geometry_arrays(model)

        self._initialize_packets(model.t_inner.value,
                                 no_of_packets)

        montecarlo.montecarlo_radial1d(
            model, plasma, self,
            virtual_packet_flag=no_of_virtual_packets,
            nthreads=nthreads,
            last_run=last_run)
        # Workaround so that j_blue_estimator is in the right ordering
        # They are written as an array of dimension (no_of_shells, no_of_lines)
        # but python expects (no_of_lines, no_of_shells)
        self.j_blue_estimator = np.ascontiguousarray(
                self.j_blue_estimator.flatten().reshape(
                    self.j_blue_estimator.shape, order='F')
                )
        self.Edotlu_estimator = np.ascontiguousarray(
                self.Edotlu_estimator.flatten().reshape(
                    self.Edotlu_estimator.shape, order='F')
                )

    def legacy_return(self):
        return (self.output_nu, self.output_energy,
                self.j_estimator, self.nu_bar_estimator,
                self.last_line_interaction_in_id,
                self.last_line_interaction_out_id,
                self.last_interaction_type,
                self.last_line_interaction_shell_id)

    def get_line_interaction_id(self, line_interaction_type):
        return ['scatter', 'downbranch', 'macroatom'].index(
            line_interaction_type)

    @property
    def output_nu(self):
        return u.Quantity(self._output_nu, u.Hz)

    @property
    def output_energy(self):
        return u.Quantity(self._output_energy, u.erg)

    @property
    def virtual_packet_nu(self):
        try:
            return u.Quantity(self.virt_packet_nus, u.Hz)
        except AttributeError:
            warnings.warn(
                    "MontecarloRunner.virtual_packet_nu:"
                    "compile with --with-vpacket-logging"
                    "to access this property", UserWarning)
            return None

    @property
    def virtual_packet_energy(self):
        try:
            return u.Quantity(self.virt_packet_energies, u.erg)
        except AttributeError:
            warnings.warn(
                    "MontecarloRunner.virtual_packet_energy:"
                    "compile with --with-vpacket-logging"
                    "to access this property", UserWarning)
            return None

    @property
    def virtual_packet_luminosity(self):
        try:
            return self.virtual_packet_energy / self.time_of_simulation
        except TypeError:
            warnings.warn(
                    "MontecarloRunner.virtual_packet_luminosity:"
                    "compile with --with-vpacket-logging"
                    "to access this property", UserWarning)
            return None

    @property
    def packet_luminosity(self):
        return self.output_energy / self.time_of_simulation

    @property
    def emitted_packet_mask(self):
        return self.output_energy >= 0

    @property
    def emitted_packet_nu(self):
        return self.output_nu[self.emitted_packet_mask]

    @property
    def reabsorbed_packet_nu(self):
        return self.output_nu[~self.emitted_packet_mask]

    @property
    def emitted_packet_luminosity(self):
        return self.packet_luminosity[self.emitted_packet_mask]

    @property
    def reabsorbed_packet_luminosity(self):
        return -self.packet_luminosity[~self.emitted_packet_mask]

    @property
    def montecarlo_reabsorbed_luminosity(self):
        return u.Quantity(
                np.histogram(
                    self.reabsorbed_packet_nu,
                    weights=self.reabsorbed_packet_luminosity,
                    bins=self.spectrum_frequency.value)[0],
                'erg / s'
                )

    @property
    def montecarlo_emitted_luminosity(self):
        return u.Quantity(
                np.histogram(
                    self.emitted_packet_nu,
                    weights=self.emitted_packet_luminosity,
                    bins=self.spectrum_frequency.value)[0],
                'erg / s'
                )

    @property
    def montecarlo_virtual_luminosity(self):
        return (
                self._montecarlo_virtual_luminosity[:-1] /
                self.time_of_simulation.value)

    def calculate_emitted_luminosity(self, luminosity_nu_start,
                                     luminosity_nu_end):

        luminosity_wavelength_filter = (
            (self.emitted_packet_nu > luminosity_nu_start) &
            (self.emitted_packet_nu < luminosity_nu_end))

        emitted_luminosity = self.emitted_packet_luminosity[
            luminosity_wavelength_filter].sum()

        return emitted_luminosity

    def calculate_reabsorbed_luminosity(
            self, luminosity_nu_start,
            luminosity_nu_end):

        luminosity_wavelength_filter = (
            (self.reabsorbed_packet_nu > luminosity_nu_start) &
            (self.reabsorbed_packet_nu < luminosity_nu_end))

        reabsorbed_luminosity = self.reabsorbed_packet_luminosity[
            luminosity_wavelength_filter].sum()

        return reabsorbed_luminosity

    def calculate_radiationfield_properties(self):
        """
        Calculate an updated radiation field from the :math:
        `\\bar{nu}_\\textrm{estimator}` and :math:`\\J_\\textrm{estimator}`
        calculated in the montecarlo simulation.
        The details of the calculation can be found in the documentation.

        Parameters
        ----------

        nubar_estimator : ~np.ndarray (float)

        j_estimator : ~np.ndarray (float)

        Returns
        -------

        t_rad : ~astropy.units.Quantity (float)

        w : ~numpy.ndarray (float)

        """

        t_rad = (
                self.t_rad_estimator_constant *
                self.nu_bar_estimator /
                self.j_estimator)
        w = self.j_estimator / (
                4 * const.sigma_sb.cgs.value * t_rad ** 4 *
                self.time_of_simulation.value *
                self.volume.value)

        return t_rad * u.K, w

    def calculate_luminosity_inner(self, model):
        return (4 * np.pi * const.sigma_sb.cgs *
                model.r_inner[0] ** 2 * model.t_inner ** 4).to('erg/s')

    def calculate_time_of_simulation(self, model):
        return (1.0 * u.erg / self.calculate_luminosity_inner(model))

    def calculate_f_nu(self, frequency):
        pass

    def calculate_f_lambda(self, wavelength):
        pass

    @classmethod
    def from_config(cls, config):
        """
        Create a new MontecarloRunner instance from a Configuration object.

        Parameters
        ----------
        config : tardis.io.config_reader.Configuration

        Returns
        -------
        MontecarloRunner

        """
        if config.plasma.disable_electron_scattering:
            logger.warn('Disabling electron scattering - this is not physical')
            sigma_thomson = 1e-200 * (u.cm ** 2)
        else:
            logger.debug("Electron scattering switched on")
            sigma_thomson = const.sigma_T.cgs

        spectrum_frequency = quantity_linspace(
            config.spectrum.stop.to('Hz', u.spectral()),
            config.spectrum.start.to('Hz', u.spectral()),
            num=config.spectrum.num + 1)

        return cls(seed=config.montecarlo.seed,
                   spectrum_frequency=spectrum_frequency,
                   virtual_spectrum_range=config.montecarlo.virtual_spectrum_range,
                   sigma_thomson=sigma_thomson,
                   enable_reflective_inner_boundary=config.montecarlo.enable_reflective_inner_boundary,
                   inner_boundary_albedo=config.montecarlo.inner_boundary_albedo,
                   line_interaction_type=config.plasma.line_interaction_type
                   )
