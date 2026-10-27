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
    VPacketCollection,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
    RPacket,
)
from tardis.transport.montecarlo.packets.trackers.tracker_full import (
    TrackerFull,
)


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
    vpacket_collections = [
        VPacketCollection(
            source_rpacket_idx,
            spectrum_frequency_grid,
            montecarlo_configuration.VPACKET_SPAWN_START_FREQUENCY,
            montecarlo_configuration.VPACKET_SPAWN_END_FREQUENCY,
            number_of_vpackets,
            number_of_vpackets,
        )
        for source_rpacket_idx in range(no_of_packets)
    ]
    return TypedList(vpacket_collections)


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

    for v_packet_idx, spectrum_idx in enumerate(v_packets_idx):
        if (
            vpacket_collection.nus[v_packet_idx] < spectrum_frequency_grid[0]
        ) or (
            vpacket_collection.nus[v_packet_idx] > spectrum_frequency_grid[-1]
        ):
            continue
        v_packets_energy_hist[spectrum_idx] += vpacket_collection.energies[
            v_packet_idx
        ]


@njit(**njit_dict)
def calculate_virtual_packet_spectrum(
    geometry_state_numba: (
        NumbaHomologousRadial1DGeometry | NumbaRadial1DGeometry
    ),
    time_explosion: float,
    opacity_state_numba: OpacityStateNumba,
    montecarlo_configuration: MonteCarloConfiguration,
    spectrum_frequency_grid: npt.NDArray[np.float64],
    trackers: TypedList[TrackerFull],
    vpacket_collections: TypedList[VPacketCollection],
    trace_vpacket_volley_function: Callable[..., None],
) -> npt.NDArray[np.float64]:
    """Calculate a virtual-packet spectrum from real-packet histories.

    Parameters
    ----------
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
    vpacket_collections : numba.typed.List
        Unused virtual packet collection for each real packet, in the same
        order as ``trackers``, e.g. from ``create_vpacket_collections``.
        Filled and finalized in place.
    trace_vpacket_volley_function
        Mode-specific virtual-packet propagation function.

    Returns
    -------
    numpy.ndarray
        Virtual-packet energy histogram on ``spectrum_frequency_grid``.

    Raises
    ------
    ValueError
        If there is not exactly one unused collection per tracker.
    """
    no_of_packets = len(trackers)
    # The collections are indexed without bounds checks and appended to in
    # place, so a mismatched or reused list (even one finalized while empty,
    # whose arrays are truncated below its length) would silently corrupt
    # the result
    if len(vpacket_collections) != no_of_packets:
        raise ValueError(
            "vpacket_collections must have one collection per tracker."
        )
    for vpacket_collection in vpacket_collections:
        if (
            vpacket_collection.idx != 0
            or len(vpacket_collection.nus) != vpacket_collection.length
        ):
            raise ValueError("vpacket_collections must be unused.")

    for i in prange(no_of_packets):
        # prange yields an unsigned index; cast to int64 so typed-list
        # indexing is not an unsafe cast and the seed sum stays integer.
        packet_idx = np.int64(i)
        tracker = trackers[packet_idx]
        vpacket_collection = vpacket_collections[packet_idx]
        # Offset the configured seed for deterministic parallel packet sampling.
        vpacket_seed = montecarlo_configuration.MONTECARLO_SEED + packet_idx
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
                last_interaction_in_r = geometry_state_numba.r_inner[0]
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

    return v_packets_energy_hist
