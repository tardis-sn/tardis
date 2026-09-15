"""Virtual-packet spectrum post-processing."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import numpy.typing as npt
from numba import njit, prange
from numba.typed import List as TypedList

from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.model.geometry.radial1d_homologous import (
    NumbaHomologousRadial1DGeometry,
)
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.transport.montecarlo import njit_dict, njit_dict_no_parallel
from tardis.transport.montecarlo.configuration.base import (
    MonteCarloConfiguration,
)
from tardis.transport.montecarlo.packets.packet_collections import (
    PacketCollection,
    VPacketCollection,
    consolidate_vpacket_tracker,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
    RPacket,
)
from tardis.transport.montecarlo.packets.trackers.tracker_full import (
    TrackerFull,
)

VPACKET_SEED_XOR = 0x80000000


@njit(**njit_dict_no_parallel)
def create_vpacket_collections(
    no_of_packets: int,
    spectrum_frequency_grid: npt.NDArray[np.float64],
    montecarlo_configuration: MonteCarloConfiguration,
    number_of_vpackets: int,
) -> TypedList[VPacketCollection]:
    """Create per-packet virtual packet collections.

    Parameters
    ----------
    no_of_packets : int
        Number of real packets.
    spectrum_frequency_grid : numpy.ndarray
        Frequency-bin edges for the real- and virtual-packet histograms.
        n + 1 edges for n bins, with uniform spacing.
    montecarlo_configuration : MonteCarloConfiguration
        Monte Carlo transport configuration.
    number_of_vpackets : int
        Number of virtual packets spawned per real packet interaction.

    Returns
    -------
    numba.typed.List
        Per-packet virtual packet collections.
    """
    vpacket_collections = TypedList()
    for i in range(no_of_packets):
        vpacket_collections.append(
            VPacketCollection(
                i,
                spectrum_frequency_grid,
                montecarlo_configuration.VPACKET_SPAWN_START_FREQUENCY,
                montecarlo_configuration.VPACKET_SPAWN_END_FREQUENCY,
                number_of_vpackets,
                montecarlo_configuration.TEMPORARY_V_PACKET_BINS,
            )
        )
    return vpacket_collections


@njit(**njit_dict_no_parallel)
def add_vpacket_collection_to_histogram(
    v_packets_energy_hist: npt.NDArray[np.float64],
    vpacket_collection: VPacketCollection,
    spectrum_frequency_grid: npt.NDArray[np.float64],
    delta_nu: float,
) -> None:
    """Add virtual packet collection energies to a spectrum histogram.

    Parameters
    ----------
    v_packets_energy_hist : numpy.ndarray
        Energy histogram to update in place.
    vpacket_collection : VPacketCollection
        Finalized virtual packet collection.
    spectrum_frequency_grid : numpy.ndarray
        Frequency-bin edges for the virtual-packet histogram.
    delta_nu : float
        Frequency bin width.
    """
    v_packets_idx = np.floor(
        (vpacket_collection.nus - spectrum_frequency_grid[0]) / delta_nu
    ).astype(np.int64)

    for j, idx in enumerate(v_packets_idx):
        if (vpacket_collection.nus[j] < spectrum_frequency_grid[0]) or (
            vpacket_collection.nus[j] > spectrum_frequency_grid[-1]
        ):
            continue
        v_packets_energy_hist[idx] += vpacket_collection.energies[j]


@njit(**njit_dict_no_parallel)
def get_vpacket_tracker(
    vpacket_collections: TypedList[VPacketCollection],
    spectrum_frequency_grid: npt.NDArray[np.float64],
    montecarlo_configuration: MonteCarloConfiguration,
) -> VPacketCollection:
    """Consolidate virtual packet tracking output when configured.

    Parameters
    ----------
    vpacket_collections : numba.typed.List
        Per-packet virtual packet collections.
    spectrum_frequency_grid : numpy.ndarray
        Frequency-bin edges for the virtual-packet histogram.
    montecarlo_configuration : MonteCarloConfiguration
        Monte Carlo transport configuration.

    Returns
    -------
    VPacketCollection
        Consolidated tracker or an empty placeholder collection.
    """
    if montecarlo_configuration.ENABLE_VPACKET_TRACKING:
        return consolidate_vpacket_tracker(
            vpacket_collections,
            spectrum_frequency_grid,
            montecarlo_configuration.VPACKET_SPAWN_START_FREQUENCY,
            montecarlo_configuration.VPACKET_SPAWN_END_FREQUENCY,
        )
    return VPacketCollection(
        -1,
        spectrum_frequency_grid,
        montecarlo_configuration.VPACKET_SPAWN_START_FREQUENCY,
        montecarlo_configuration.VPACKET_SPAWN_END_FREQUENCY,
        -1,
        1,
    )


@njit(**njit_dict)
def calculate_virtual_packet_spectrum(
    packet_collection: PacketCollection,
    geometry_state_numba: (
        NumbaHomologousRadial1DGeometry | NumbaRadial1DGeometry
    ),
    time_explosion: float,
    opacity_state_numba: OpacityStateNumba,
    montecarlo_configuration: MonteCarloConfiguration,
    spectrum_frequency_grid: npt.NDArray[np.float64],
    trackers: TypedList[TrackerFull],
    number_of_vpackets: int,
    trace_vpacket_volley_function: Callable[..., None],
) -> tuple[npt.NDArray[np.float64], VPacketCollection]:
    """Calculate a virtual-packet spectrum from real-packet histories.

    Parameters
    ----------
    packet_collection : PacketCollection
        Real-packet inputs, including the per-packet random seeds.
    geometry_state_numba
        Numba geometry object used for real-packet transport.
    time_explosion : float
        Time since explosion in seconds. Non-homologous mode accepts but does
        not use this value.
    opacity_state_numba
        Frozen opacity state used for real-packet transport.
    montecarlo_configuration : MonteCarloConfiguration
        Monte Carlo transport configuration.
    spectrum_frequency_grid : numpy.ndarray
        Frequency-bin edges for the virtual-packet histogram.
    trackers
        Completed full-event tracker for each real packet.
    number_of_vpackets : int
        Number of virtual packets spawned per launch or interaction.
    trace_vpacket_volley_function
        Mode-specific virtual-packet propagation function.

    Returns
    -------
    tuple
        Virtual-packet energy histogram and consolidated logging collection.
    """
    no_of_packets = len(trackers)
    vpacket_collections = create_vpacket_collections(
        no_of_packets,
        spectrum_frequency_grid,
        montecarlo_configuration,
        number_of_vpackets,
    )

    for i in prange(no_of_packets):
        packet_idx = np.int64(i)
        tracker = trackers[packet_idx]
        vpacket_collection = vpacket_collections[packet_idx]
        vpacket_seed = (
            packet_collection.packet_seeds[packet_idx] ^ VPACKET_SEED_XOR
        )
        np.random.seed(vpacket_seed)

        for event_idx in range(tracker.event_id):
            interaction_type = tracker.interaction_type[event_idx]
            is_initial_launch = (
                interaction_type == InteractionType.BOUNDARY
                and tracker.before_shell_id[event_idx] == -1
                and tracker.after_shell_id[event_idx] == 0
            )
            is_physical_interaction = (
                interaction_type == InteractionType.LINE
                or interaction_type == InteractionType.ESCATTERING
            )
            if not (is_initial_launch or is_physical_interaction):
                continue

            r_packet = RPacket(
                tracker.radius[event_idx],
                tracker.after_mu[event_idx],
                tracker.after_nu[event_idx],
                tracker.after_energy[event_idx],
                vpacket_seed,
                packet_idx,
            )
            r_packet.current_shell_id = tracker.after_shell_id[event_idx]
            r_packet.next_line_id = tracker.next_line_id[event_idx]

            if is_initial_launch:
                last_interaction_in_nu = 0.0
                last_interaction_in_r = 0.0
                last_interaction_type = InteractionType.NO_INTERACTION
                last_interaction_in_id = -1
                last_interaction_out_id = -1
                last_interaction_shell_id = -1
            else:
                last_interaction_in_nu = tracker.before_nu[event_idx]
                last_interaction_in_r = tracker.radius[event_idx]
                last_interaction_type = interaction_type
                last_interaction_in_id = tracker.line_absorb_id[event_idx]
                last_interaction_out_id = tracker.line_emit_id[event_idx]
                last_interaction_shell_id = tracker.before_shell_id[event_idx]

            trace_vpacket_volley_function(
                r_packet,
                vpacket_collection,
                geometry_state_numba,
                time_explosion,
                opacity_state_numba,
                montecarlo_configuration.ENABLE_FULL_RELATIVITY,
                montecarlo_configuration.VPACKET_TAU_RUSSIAN,
                montecarlo_configuration.SURVIVAL_PROBABILITY,
                last_interaction_in_nu,
                last_interaction_in_r,
                last_interaction_type,
                last_interaction_in_id,
                last_interaction_out_id,
                last_interaction_shell_id,
            )

        vpacket_collection.finalize_arrays()

    v_packets_energy_hist = np.zeros_like(spectrum_frequency_grid)
    delta_nu = spectrum_frequency_grid[1] - spectrum_frequency_grid[0]
    for vpacket_collection in vpacket_collections:
        add_vpacket_collection_to_histogram(
            v_packets_energy_hist,
            vpacket_collection,
            spectrum_frequency_grid,
            delta_nu,
        )

    vpacket_tracker = get_vpacket_tracker(
        vpacket_collections,
        spectrum_frequency_grid,
        montecarlo_configuration,
    )
    return v_packets_energy_hist, vpacket_tracker
