import os
import logging
import warnings

from astropy import units as u, constants as const

from scipy.special import zeta
from spectrum import TARDISSpectrum

from tardis.montecarlo import montecarlo, packet_source
from tardis.io.util import to_hdf

import numpy as np
import pandas as pd

logger = logging.getLevelName(__name__)

class MontecarloRunner(object):
    """
    This class is designed as an interface between the Python part and the
    montecarlo C-part
    """

    w_estimator_constant = ((const.c ** 2 / (2 * const.h)) *
                            (15 / np.pi ** 4) * (const.h / const.k_B) ** 4 /
                            (4 * np.pi)).cgs.value

    t_rad_estimator_constant = ((np.pi**4 / (15 * 24 * zeta(5, 1))) *
                                (const.h / const.k_B)).cgs.value

    def __init__(self, seed, spectrum_frequency, distance=None):
        self.packet_source = packet_source.BlackBodySimpleSource(seed)
        self.spectrum_frequency = spectrum_frequency
        self.spectrum = TARDISSpectrum(spectrum_frequency, distance)
        self.spectrum_virtual = TARDISSpectrum(spectrum_frequency, distance)
        self.spectrum_reabsorbed = TARDISSpectrum(spectrum_frequency, distance)


    def _initialize_estimator_arrays(self, no_of_shells, tau_sobolev_shape):
        """
        Initialize the output arrays of the montecarlo simulation.

        Parameters
        ----------

        model: ~Radial1DModel
        """

        #Estimators
        self.j_estimator = np.zeros(no_of_shells, dtype=np.float64)
        self.nu_bar_estimator = np.zeros(no_of_shells, dtype=np.float64)
        self.j_blue_estimator = np.zeros(tau_sobolev_shape)
        self.Edotlu_estimator = np.zeros(tau_sobolev_shape)


    def _initialize_geometry_arrays(self, structure):
        """
        Generate the cgs like geometry arrays for the montecarlo part

        Parameters
        ----------

        structure: ~ConfigurationNameSpace
        """
        self.r_inner_cgs = structure.r_inner.to('cm').value
        self.r_outer_cgs = structure.r_outer.to('cm').value
        self.v_inner_cgs = structure.v_inner.to('cm/s').value

    def _initialize_packets(self, T, no_of_packets, no_of_virtual_packets=None):
        nus, mus, energies = self.packet_source.create_packets(T, no_of_packets)
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
        self.last_interaction_type = -1 * np.ones(no_of_packets, dtype=np.int64)
        self.last_interaction_in_nu = np.zeros(no_of_packets, dtype=np.float64)


        self.legacy_montecarlo_virtual_luminosity = np.zeros_like(
            self.spectrum_frequency.value)

    def legacy_update_spectrum(self, no_of_virtual_packets):
        montecarlo_reabsorbed_luminosity = np.histogram(
            self.reabsorbed_packet_nu,
            weights=self.reabsorbed_packet_luminosity,
            bins=self.spectrum_frequency.value)[0] * u.erg / u.s

        montecarlo_emitted_luminosity = np.histogram(
            self.emitted_packet_nu,
            weights=self.emitted_packet_luminosity,
            bins=self.spectrum_frequency.value)[0] * u.erg / u.s

        self.spectrum.update_luminosity(montecarlo_emitted_luminosity)
        self.spectrum_reabsorbed.update_luminosity(
            montecarlo_reabsorbed_luminosity)

        if no_of_virtual_packets > 0:
            self.montecarlo_virtual_luminosity = (
                self.legacy_montecarlo_virtual_luminosity *
                1 * u.erg / self.time_of_simulation)[:-1]
            self.spectrum_virtual.update_luminosity(
                self.montecarlo_virtual_luminosity)

    def run(self, model, no_of_packets, no_of_virtual_packets=0, nthreads=1,last_run=False):
        """
        Running the TARDIS simulation

        Parameters
        ----------

        :param model:
        :param no_of_virtual_packets:
        :param nthreads:
        :return:
        """
        self.time_of_simulation = model.time_of_simulation
        self.volume = model.tardis_config.structure.volumes
        self._initialize_estimator_arrays(self.volume.shape[0],
                                          model.plasma.tau_sobolevs.shape)
        self._initialize_geometry_arrays(model.tardis_config.structure)

        self._initialize_packets(model.t_inner.value,
                                 no_of_packets)

        montecarlo.montecarlo_radial1d(
            model, self, virtual_packet_flag=no_of_virtual_packets,
            nthreads=nthreads,last_run=last_run)
        # Workaround so that j_blue_estimator is in the right ordering
        # They are written as an array of dimension (no_of_shells, no_of_lines)
        # but python expects (no_of_lines, no_of_shells)
        self.j_blue_estimator = np.ascontiguousarray(
                self.j_blue_estimator.flatten().reshape(
                self.j_blue_estimator.shape, order='F')
                )
        self.Edotlu_estimator = self.Edotlu_estimator.flatten().reshape(
                self.Edotlu_estimator.shape, order='F')

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
            warnings.warn("MontecarloRunner.virtual_packet_nu:"
                    "compile with --with-vpacket-logging"
                    "to access this property", UserWarning)
            return None

    @property
    def virtual_packet_energy(self):
        try:
            return u.Quantity(self.virt_packet_energies, u.erg)
        except AttributeError:
            warnings.warn("MontecarloRunner.virtual_packet_energy:"
                    "compile with --with-vpacket-logging"
                    "to access this property", UserWarning)
            return None

    @property
    def virtual_packet_luminosity(self):
        try:
            return self.virtual_packet_energy / self.time_of_simulation
        except TypeError:
            warnings.warn("MontecarloRunner.virtual_packet_luminosity:"
                    "compile with --with-vpacket-logging"
                    "to access this property", UserWarning)
            return None

    @property
    def packet_luminosity(self):
        return self.output_energy / self.time_of_simulation

    @property
    def emitted_packet_mask(self):
        return self.output_energy >=0

    @property
    def emitted_packet_nu(self):
        return self.output_nu[self.emitted_packet_mask]

    @property
    def reabsorbed_packet_nu(self):
        return self.output_nu[~self.emitted_packet_mask]

    @property
    def reabsorbed_packet_luminosity(self):
        return -self.packet_luminosity[~self.emitted_packet_mask]


    @property
    def emitted_packet_luminosity(self):
        return self.packet_luminosity[self.emitted_packet_mask]

    @property
    def reabsorbed_packet_luminosity(self):
        return -self.packet_luminosity[~self.emitted_packet_mask]

    def calculate_emitted_luminosity(self, luminosity_nu_start,
                                     luminosity_nu_end):

        luminosity_wavelength_filter = (
            (self.emitted_packet_nu > luminosity_nu_start) &
            (self.emitted_packet_nu < luminosity_nu_end))

        emitted_luminosity = self.emitted_packet_luminosity[
            luminosity_wavelength_filter].sum()

        return emitted_luminosity

    def calculate_reabsorbed_luminosity(self, luminosity_nu_start,
                                     luminosity_nu_end):

        luminosity_wavelength_filter = (
            (self.reabsorbed_packet_nu > luminosity_nu_start) &
            (self.reabsorbed_packet_nu < luminosity_nu_end))

        reabsorbed_luminosity = self.reabsorbed_packet_luminosity[
            luminosity_wavelength_filter].sum()

        return reabsorbed_luminosity


    def calculate_radiationfield_properties(self):
        """
        Calculate an updated radiation field from the :math:`\\bar{nu}_\\textrm{estimator}` and :math:`\\J_\\textrm{estimator}`
        calculated in the montecarlo simulation. The details of the calculation can be found in the documentation.

        Parameters
        ----------

        nubar_estimator : ~np.ndarray (float)

        j_estimator : ~np.ndarray (float)

        Returns
        -------

        t_rad : ~astropy.units.Quantity (float)

        w : ~numpy.ndarray (float)

        """


        t_rad = (self.t_rad_estimator_constant * self.nu_bar_estimator
                / self.j_estimator)
        w = self.j_estimator / (4 * const.sigma_sb.cgs.value * t_rad ** 4
                                * self.time_of_simulation.value
                                * self.volume.value)

        return t_rad * u.K, w

    def calculate_f_nu(self, frequency):
        pass

    def calculate_f_lambda(self, wavelength):
        pass

    def to_hdf(self, path_or_buf, path=''):
        """
        Store the runner to an HDF structure.

        Parameters
        ----------
        path_or_buf:
            Path or buffer to the HDF store
        path:
            Path inside the HDF store to store the runner
        Returns
        -------
            : None

        """
        runner_path = os.path.join(path, 'runner')
        properties = ['output_nu', 'output_energy', 'nu_bar_estimator',
                      'j_estimator', 'montecarlo_virtual_luminosity']
        to_hdf(path_or_buf, runner_path, {name: getattr(self, name) for name
                                          in properties})
        self.spectrum.to_hdf(path_or_buf, runner_path)
        self.spectrum_virtual.to_hdf(path_or_buf, runner_path,
                                     'spectrum_virtual')
        self.spectrum_reabsorbed.to_hdf(path_or_buf, runner_path,
                                        'spectrum_reabsorbed')
