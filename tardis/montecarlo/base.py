import os
import logging
import warnings
import h5py
from astropy import units as u, constants as const

from scipy.special import zeta
from spectrum import TARDISSpectrum

from tardis.util import quantity_linspace
from tardis.montecarlo import montecarlo, packet_source
from tardis.io.util import to_hdf
from tardis.io.config_reader import ConfigurationNameSpace
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

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

    def __init__(self, seed, spectrum_frequency, virtual_spectrum_range,
                 sigma_thomson, enable_reflective_inner_boundary,
                 inner_boundary_albedo, line_interaction_type, distance=None):

        self.seed = seed
        self.packet_source = packet_source.BlackBodySimpleSource(seed)
        self.spectrum_frequency = spectrum_frequency
        self.virtual_spectrum_range = virtual_spectrum_range
        self.sigma_thomson = sigma_thomson
        self.enable_reflective_inner_boundary = enable_reflective_inner_boundary
        self.inner_boundary_albedo = inner_boundary_albedo
        self.line_interaction_type = line_interaction_type
        self.distance = distance
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

    def run(self, model, plasma, no_of_packets, no_of_virtual_packets=0, nthreads=1,last_run=False):
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
        self.time_of_simulation = self.calculate_time_of_simulation(model)
        self.volume = model.volume
        self._initialize_estimator_arrays(self.volume.shape[0],
                                          plasma.tau_sobolevs.shape)
        self._initialize_geometry_arrays(model)

        self._initialize_packets(model.t_inner.value,
                                 no_of_packets)

        montecarlo.montecarlo_radial1d(
            model, plasma, self, virtual_packet_flag=no_of_virtual_packets,
            nthreads=nthreads,last_run=last_run)
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

    def calculate_luminosity_inner(self, model):
        return (4 * np.pi * const.sigma_sb.cgs *
                model.r_inner[0] ** 2 * model.t_inner ** 4).to('erg/s')

    def calculate_time_of_simulation(self, model):
        return (1.0 * u.erg / self.calculate_luminosity_inner(model))

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
                      'j_estimator', 'montecarlo_virtual_luminosity',
                      'last_interaction_in_nu',
                      'last_line_interaction_in_id',
                      'last_line_interaction_out_id',
                      'last_line_interaction_shell_id',
                      'seed', 'spectrum_frequency',
                      'virtual_spectrum_range', 'sigma_thomson', 'enable_reflective_inner_boundary',
                      'inner_boundary_albedo', 'line_interaction_type', 'distance'
                      ]
        to_hdf(path_or_buf, runner_path, {name: getattr(self, name) for name
                                          in properties})
        self.spectrum.to_hdf(path_or_buf, runner_path)
        self.spectrum_virtual.to_hdf(path_or_buf, runner_path,
                                     'spectrum_virtual')
        self.spectrum_reabsorbed.to_hdf(path_or_buf, runner_path,
                                        'spectrum_reabsorbed')

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
            sigma_thomson = 1e-200 / (u.cm ** 2)
        else:
            logger.debug("Electron scattering switched on")
            sigma_thomson = 6.652486e-25 / (u.cm ** 2)

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
                   line_interaction_type=config.plasma.line_interaction_type,
                   distance=config.supernova.get('distance', None))

    @classmethod
    def from_hdf(cls, path, h5_file, file_path,model,plasma):
        """
        This function returns a MontecarloRunner object 
        from given HDF5 File.

        Parameters
        ----------
        path : 'str'
            Path to transverse in hdf file
        h5_file : 'h5py.File'
            Given HDF5 file
        file_path : 'str'
            Path of Simulation generated HDF file 

        Returns
        -------
        `~MontecarloRunner`
        """

        if not h5_file:
            raise ValueError("h5_file Parameter can`t be None")

        runner_path = path + '/runner'
        runner_keys = ['scalars', 'spectrum_frequency',
                       'virtual_spectrum_range', 'distance','packet_luminosity','output_energy','output_nu',
                       'last_line_interaction_in_id','last_interaction_in_nu','last_line_interaction_out_id','last_line_interaction_shell_id',
                       'j_estimator','montecarlo_virtual_luminosity','nu_bar_estimator']
        runner_dict = {}
        with pd.HDFStore(file_path, 'r') as data:
            for key in h5_file[runner_path].keys():
                if key in runner_keys:
                    runner_dict[key] = {}
                    buff_path = runner_path + '/' + key + '/'
                    runner_dict[key] = data[buff_path]

        #Creates corresponding astropy.units.Quantity objects

        seed = runner_dict['scalars']['seed']
        sigma_thomson = runner_dict['scalars']['sigma_thomson'] * (1/(u.cm * u.cm))
        enable_reflective_inner_boundary = runner_dict['scalars']['enable_reflective_inner_boundary']
        inner_boundary_albedo = runner_dict['scalars']['inner_boundary_albedo']
        line_interaction_type = runner_dict['scalars']['line_interaction_type']
        distance = runner_dict['distance']
        spectrum_frequency = np.array(runner_dict['spectrum_frequency']) * u.Hz
        virtual_spectrum_range = dict(
            stop=runner_dict['virtual_spectrum_range']['stop'][0],
            start=runner_dict['virtual_spectrum_range']['start'][0],
            num=runner_dict['virtual_spectrum_range']['num'][0])
        virtual_spectrum_range = ConfigurationNameSpace(virtual_spectrum_range)

        runner =  cls(seed, spectrum_frequency, virtual_spectrum_range,
                   sigma_thomson, enable_reflective_inner_boundary,
                   inner_boundary_albedo, line_interaction_type, distance)
        
        runner.time_of_simulation = runner.calculate_time_of_simulation(model)
        runner.volume = model.volume
        runner._initialize_estimator_arrays(runner.volume.shape[0],
                                          plasma.tau_sobolevs.shape)
        runner._initialize_geometry_arrays(model)

        consts_path = path + '/consts'
        consts_keys = ['scalars']
        consts = {}
        with pd.HDFStore(file_path, 'r') as data:
            for key in h5_file[consts_path].keys():
                if key in consts_keys:
                    consts[key] = {}
                    buff_path = consts_path + '/' + key + '/'
                    consts[key] = data[buff_path]
        
        runner._initialize_packets(model.t_inner.value,
                                 int(consts['scalars']['no_of_packets']))
                                 
        # runner._initialize_packets(model.t_inner.value,
        #                          100000)
        #print int(consts['scalars']['no_of_packets'])
        # runner.j_blue_estimator = np.ascontiguousarray(
        #         runner.j_blue_estimator.flatten().reshape(
        #         runner.j_blue_estimator.shape, order='F')
        #         )
        # runner.Edotlu_estimator = np.ascontiguousarray(
        #         runner.Edotlu_estimator.flatten().reshape(
        #         runner.Edotlu_estimator.shape, order='F')
        #         )
        # runner.virt_logging = 0
        # montecarlo.montecarlo_radial1d(
        #     model, plasma, runner, virtual_packet_flag=consts['scalars']['no_of_virtual_packets'],nthreads=consts['scalars']['nthreads'],last_run=True)
        runner._output_energy = runner_dict['output_energy']
        runner._output_nu = runner_dict['output_nu']
        runner.last_line_interaction_in_id = np.array(runner_dict['last_line_interaction_in_id'])
        runner.last_interaction_in_nu = np.array(runner_dict['last_interaction_in_nu'])
        runner.last_line_interaction_out_id = np.array(runner_dict['last_line_interaction_out_id'])
        runner.last_line_interaction_shell_id = np.array(runner_dict['last_line_interaction_shell_id'])
        runner.j_estimator = np.array(runner_dict['j_estimator'])
        runner.montecarlo_virtual_luminosity = np.array(runner_dict['montecarlo_virtual_luminosity']) * u.erg/u.s
        runner.nu_bar_estimator = np.array(runner_dict['nu_bar_estimator'])

        runner.line_lists_tau_sobolevs =  plasma.tau_sobolevs.values.flatten(order='F')
        if runner.get_line_interaction_id(runner.line_interaction_type)>=1:
            runner.transition_probabilities = (
                  plasma.transition_probabilities.values.flatten(order='F'))
        
        # runner.virt_packet_nus = np.zeros(0)
        # runner.virt_packet_energies = np.zeros(0)             
        # runner.virt_packet_last_interaction_in_nu = np.zeros(0)
        # runner.virt_packet_last_interaction_type = np.zeros(0)
        runner.inverse_electron_densities = (1.0 / plasma.electron_densities.values)
        #print runner.get_line_interaction_id(runner.line_interaction_type)
        
        #see no_of_packets
        #print runner.sigma_thomson

        #runner.run(model, plasma, int(consts['scalars']['no_of_packets']), no_of_virtual_packets=int(consts['scalars']['no_of_virtual_packets']), nthreads=int(consts['scalars']['nthreads']),last_run=False)

        return runner