import warnings

import numpy as np
from astropy import units as u

from tardis.io.util import HDFWriterMixin
from tardis.transport.montecarlo.estimators.dilute_blackbody_properties import (
    MCDiluteBlackBodyRadFieldSolver,
)
from tardis.spectrum.formal_integral import IntegrationError
from tardis.spectrum.spectrum import TARDISSpectrum


class MonteCarloTransportState(HDFWriterMixin):
    hdf_properties = [
        "output_nu",
        "output_energy",
        "nu_bar_estimator",
        "j_estimator",
        "j_blue_estimator",
        "packet_luminosity",
        "time_of_simulation",
        "emitted_packet_mask",
        "last_interaction_type",
        "last_interaction_in_nu",
        "last_interaction_in_r",
        "last_line_interaction_out_id",
        "last_line_interaction_in_id",
        "last_line_interaction_shell_id",
    ]

    vpacket_hdf_properties = [
        "virt_packet_nus",
        "virt_packet_energies",
        "virt_packet_initial_rs",
        "virt_packet_initial_mus",
        "virt_packet_last_interaction_in_nu",
        "virt_packet_last_interaction_in_r",
        "virt_packet_last_interaction_type",
        "virt_packet_last_line_interaction_in_id",
        "virt_packet_last_line_interaction_out_id",
        "virt_packet_last_line_interaction_shell_id",
    ]

    hdf_name = "transport_state"

    last_interaction_type = None
    last_interaction_in_nu = None
    last_interaction_in_r = None
    last_line_interaction_out_id = None
    last_line_interaction_in_id = None
    last_line_interaction_shell_id = None

    virt_logging = False

    def __init__(
        self,
        packet_collection,
        radfield_mc_estimators,
        geometry_state,
        opacity_state,
        rpacket_tracker=None,
        vpacket_tracker=None,
    ):
        self.packet_collection = packet_collection
        self.radfield_mc_estimators = radfield_mc_estimators
        self.enable_full_relativity = False
        self.enable_continuum_processes = False
        self.geometry_state = geometry_state
        self.opacity_state = opacity_state
        self.rpacket_tracker = rpacket_tracker
        self.vpacket_tracker = vpacket_tracker

    def calculate_radiationfield_properties(self):
        """
        Calculate an updated radiation field from the :math:
        `\\bar{nu}_\\textrm{estimator}` and :math:`\\J_\\textrm{estimator}`
        calculated in the montecarlo simulation.
        The details of the calculation can be found in the documentation.

        Parameters
        ----------
        nubar_estimator : np.ndarray (float)
        j_estimator : np.ndarray (float)

        Returns
        -------
        t_radiative : astropy.units.Quantity (float)
        dilution_factor : numpy.ndarray (float)
        """
        dilute_bb_solver = MCDiluteBlackBodyRadFieldSolver()
        dilute_bb_radfield = dilute_bb_solver.solve(
            self.radfield_mc_estimators,
            self.time_of_simulation,
            self.geometry_state.volume,
        )

        return (
            dilute_bb_radfield.temperature,
            dilute_bb_radfield.dilution_factor,
        )

    @property
    def output_nu(self):
        return self.packet_collection.output_nus * u.Hz

    @property
    def output_energy(self):
        return self.packet_collection.output_energies * u.erg

    @property
    def nu_bar_estimator(self):
        return self.radfield_mc_estimators.nu_bar_estimator

    @property
    def j_estimator(self):
        return self.radfield_mc_estimators.j_estimator

    @property
    def j_blue_estimator(self):
        return self.radfield_mc_estimators.j_blue_estimator

    @property
    def time_of_simulation(self):
        return self.packet_collection.time_of_simulation * u.s

    @property
    def packet_luminosity(self):
        return (
            self.packet_collection.output_energies
            * u.erg
            / (self.packet_collection.time_of_simulation * u.s)
        )

    @property
    def emitted_packet_mask(self):
        return self.packet_collection.output_energies >= 0

    @property
    def emitted_packet_nu(self):
        return (
            self.packet_collection.output_nus[self.emitted_packet_mask] * u.Hz
        )

    @property
    def reabsorbed_packet_nu(self):
        return (
            self.packet_collection.output_nus[~self.emitted_packet_mask] * u.Hz
        )

    @property
    def emitted_packet_luminosity(self):
        return self.packet_luminosity[self.emitted_packet_mask]

    @property
    def reabsorbed_packet_luminosity(self):
        return -self.packet_luminosity[~self.emitted_packet_mask]

    @property
    def virt_packet_nus(self):
        try:
            return u.Quantity(self.vpacket_tracker.nus, u.Hz)
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virt_packet_nus:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virt_packet_energies(self):
        try:
            return u.Quantity(self.vpacket_tracker.energies, u.erg)
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virt_packet_energies:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virtual_packet_luminosity(self):
        try:
            return (
                self.virt_packet_energies
                / self.packet_collection.time_of_simulation
            )
        except TypeError:
            warnings.warn(
                "MontecarloTransport.virtual_packet_luminosity:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virt_packet_initial_rs(self):
        try:
            return u.Quantity(self.vpacket_tracker.initial_rs, u.erg)
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virt_packet_initial_rs:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virt_packet_initial_mus(self):
        try:
            return u.Quantity(self.vpacket_tracker.initial_mus, u.erg)
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virt_packet_initial_mus:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virt_packet_last_interaction_in_nu(self):
        try:
            return u.Quantity(self.vpacket_tracker.last_interaction_in_nu, u.Hz)
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virt_packet_last_interaction_in_nu:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virt_packet_last_interaction_in_r(self):
        try:
            return u.Quantity(self.vpacket_tracker.last_interaction_in_r, u.cm)
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virt_packet_last_interaction_in_r:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virt_packet_last_interaction_type(self):
        try:
            return self.vpacket_tracker.last_interaction_type
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virt_packet_last_interaction_type:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virt_packet_last_line_interaction_in_id(self):
        try:
            return self.vpacket_tracker.last_interaction_in_id
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virt_packet_last_line_interaction_in_id:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virt_packet_last_line_interaction_out_id(self):
        try:
            return self.vpacket_tracker.last_interaction_out_id
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virt_packet_last_line_interaction_out_id:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virt_packet_last_line_interaction_shell_id(self):
        try:
            return self.vpacket_tracker.last_interaction_shell_id
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virt_packet_last_line_interaction_shell_id:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None
