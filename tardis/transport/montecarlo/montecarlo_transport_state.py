import warnings

import pandas as pd
from astropy import units as u

from tardis.io.hdf_writer_mixin import HDFWriterMixin
from tardis.model.geometry.radial1d import (
    HomologousRadial1DGeometry,
    NumbaRadial1DGeometry,
)
from tardis.opacities.macro_atom.macroatom_state import MacroAtomState
from tardis.opacities.opacity_state import OpacityState
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.plasma.base import BasePlasma
from tardis.transport.montecarlo.estimators.radfield_mc_estimators import (
    RadiationFieldMCEstimators,
    initialize_estimator_statistics,
)
from tardis.transport.montecarlo.packet_source.base import BasePacketSource
from tardis.transport.montecarlo.packets.packet_collections import (
    PacketCollection,
    VPacketCollection,
)


class MonteCarloTransportState(HDFWriterMixin):
    """
    Monte Carlo transport state containing packet data and radiation field estimators.

    Attributes
    ----------
    packet_collection : PacketCollection
        Collection of Monte Carlo packets with initial and final properties.
    radfield_mc_estimators : RadiationFieldMCEstimators
        Estimators for radiation field properties.
    geometry_state : NumbaRadial1DGeometry
        Numba-compiled geometry state with shell boundaries.
    opacity_state : OpacityStateNumba
        Numba-compiled opacity state with line and continuum opacities.
    time_explosion : u.Quantity
        Time since explosion.
    enable_full_relativity : bool
        Flag for full relativistic corrections.
    enable_continuum_processes : bool
        Flag for continuum processes.
    tracker_full_df : pd.DataFrame | None
        Full packet tracking data.
    tracker_last_interaction_df : pd.DataFrame | None
        Last interaction tracking data for each packet.
    vpacket_tracker : VPacketCollection | None
        Virtual packet tracking collection.
    last_interaction_type : None
        Last interaction type placeholder.
    last_interaction_in_nu : None
        Last interaction frequency placeholder.
    last_interaction_in_r : None
        Last interaction radius placeholder.
    last_line_interaction_out_id : None
        Last line interaction output ID placeholder.
    last_line_interaction_in_id : None
        Last line interaction input ID placeholder.
    last_line_interaction_shell_id : None
        Last line interaction shell ID placeholder.
    virt_logging : bool
        Flag for virtual packet logging.
    """

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

    hdf_name: str = "transport_state"

    packet_collection: PacketCollection
    radfield_mc_estimators: RadiationFieldMCEstimators
    geometry_state: NumbaRadial1DGeometry
    opacity_state: OpacityStateNumba
    time_explosion: u.Quantity
    enable_full_relativity: bool
    enable_continuum_processes: bool
    tracker_full_df: pd.DataFrame | None
    tracker_last_interaction_df: pd.DataFrame | None
    vpacket_tracker: VPacketCollection | None

    last_interaction_type = None
    last_interaction_in_nu = None
    last_interaction_in_r = None
    last_line_interaction_out_id = None
    last_line_interaction_in_id = None
    last_line_interaction_shell_id = None

    virt_logging: bool = False

    def __init__(
        self,
        packet_collection: PacketCollection,
        radfield_mc_estimators: RadiationFieldMCEstimators,
        geometry_state: NumbaRadial1DGeometry,
        opacity_state: OpacityStateNumba,
        time_explosion: u.Quantity,
        tracker_full_df: pd.DataFrame | None = None,
        tracker_last_interaction_df: pd.DataFrame | None = None,
        vpacket_tracker: VPacketCollection | None = None,
    ) -> None:
        """
        Initialize a Monte Carlo transport state.

        Parameters
        ----------
        packet_collection : PacketCollection
            Collection of Monte Carlo packets with initial properties.
        radfield_mc_estimators : RadiationFieldMCEstimators
            Estimators for radiation field properties.
        geometry_state : NumbaRadial1DGeometry
            Numba-compiled geometry state with shell boundaries.
        opacity_state : OpacityStateNumba
            Numba-compiled opacity state with line and continuum opacities.
        time_explosion : u.Quantity
            Time since explosion.
        tracker_full_df : pd.DataFrame, optional
            Full packet tracking data. Default is None.
        tracker_last_interaction_df : pd.DataFrame, optional
            Last interaction tracking data for each packet. Default is None.
        vpacket_tracker : VPacketCollection, optional
            Virtual packet tracking collection. Default is None.
        """
        self.packet_collection = packet_collection
        self.radfield_mc_estimators = radfield_mc_estimators
        self.enable_full_relativity = False
        self.enable_continuum_processes = False
        self.time_explosion = time_explosion
        self.geometry_state = geometry_state
        self.opacity_state = opacity_state
        self.tracker_full_df = tracker_full_df
        self.tracker_last_interaction_df = tracker_last_interaction_df
        self.vpacket_tracker = vpacket_tracker

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


def initialize_transport_state(
    geometry_state: HomologousRadial1DGeometry,
    time_explosion: u.Quantity,
    opacity_state: OpacityState,
    macro_atom_state: MacroAtomState,
    plasma: BasePlasma,
    packet_source: BasePacketSource,
    no_of_packets: int,
    line_interaction_type: str,
    iteration: int = 0,
) -> MonteCarloTransportState:
    """
    Initialize a MonteCarloTransportState from simulation components.

    Parameters
    ----------
    geometry_state : HomologousRadial1DGeometry
        Geometry state containing shell boundaries and velocity information.
    time_explosion : u.Quantity
        Time since explosion.
    opacity_state : OpacityState
        Opacity state containing line and continuum opacities.
    macro_atom_state : MacroAtomState
        Macro atom state with transition probabilities.
    plasma : BasePlasma
        Plasma state with physical properties.
    packet_source : BasePacketSource
        Packet source for creating initial packets.
    no_of_packets : int
        Number of packets to create.
    line_interaction_type : str
        Type of line interaction (e.g., 'scatter', 'macroatom', 'downbranch').
    iteration : int, optional
        Iteration number for seed offset. Default is 0.

    Returns
    -------
    MonteCarloTransportState
        Initialized transport state ready for Monte Carlo simulation.
    """
    if not plasma.continuum_interaction_species.empty:
        gamma_shape = plasma.gamma.shape
    else:
        gamma_shape = (0, 0)

    packet_collection = packet_source.create_packets(
        no_of_packets, seed_offset=iteration
    )

    geometry_state_numba = geometry_state.to_numba()
    opacity_state_numba = opacity_state.to_numba(
        macro_atom_state,
        line_interaction_type,
    )
    opacity_state_numba = opacity_state_numba[
        geometry_state.v_inner_boundary_index : geometry_state.v_outer_boundary_index
    ]

    estimators = initialize_estimator_statistics(
        opacity_state_numba.tau_sobolev.shape, gamma_shape
    )

    return MonteCarloTransportState(
        packet_collection,
        estimators,
        geometry_state=geometry_state_numba,
        opacity_state=opacity_state_numba,
        time_explosion=time_explosion,
    )
