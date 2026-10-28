"""Virtual-packet spectrum post-processing."""

from __future__ import annotations

from collections.abc import Callable
from typing import NamedTuple

import numpy as np
import numpy.typing as npt
import pandas as pd
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
from tardis.transport.montecarlo.packets.trackers import (
    interaction_type_to_codes,
)

TRACKER_FULL_DF_COLUMNS = [
    "interaction_type",
    "radius",
    "before_shell_id",
    "after_shell_id",
    "before_nu",
    "after_nu",
    "after_mu",
    "after_energy",
    "line_absorb_id",
    "line_emit_id",
    "next_line_id",
]


class TrackerEventArrays(NamedTuple):
    """
    Real-packet event histories as flat arrays for virtual-packet replay.

    Events of packet ``packet_id[i]`` occupy
    ``event_offsets[i]:event_offsets[i + 1]`` in every per-event array, in
    event order.
    """

    packet_id: npt.NDArray[np.int64]
    event_offsets: npt.NDArray[np.int64]
    interaction_type: npt.NDArray[np.int64]
    radius: npt.NDArray[np.float64]
    before_shell_id: npt.NDArray[np.int64]
    after_shell_id: npt.NDArray[np.int64]
    before_nu: npt.NDArray[np.float64]
    after_nu: npt.NDArray[np.float64]
    after_mu: npt.NDArray[np.float64]
    after_energy: npt.NDArray[np.float64]
    line_absorb_id: npt.NDArray[np.int64]
    line_emit_id: npt.NDArray[np.int64]
    next_line_id: npt.NDArray[np.int64]


def tracker_full_df_to_event_arrays(
    tracker_full_df: pd.DataFrame,
) -> TrackerEventArrays:
    """
    Extract the event arrays needed to replay real-packet histories.

    Parameters
    ----------
    tracker_full_df : pandas.DataFrame
        Full real-packet tracking DataFrame indexed by ``packet_id`` and
        ``event_id``, e.g. ``transport_state.tracker_full_df`` or one loaded
        from disk.

    Returns
    -------
    TrackerEventArrays
        Flat event arrays grouped by packet.

    Raises
    ------
    ValueError
        If the ``packet_id`` index level or a column needed to replay the
        histories is missing.
    """
    if "packet_id" not in tracker_full_df.index.names:
        raise ValueError(
            "tracker_full_df must be indexed by packet_id and event_id."
        )
    missing_columns = [
        column
        for column in TRACKER_FULL_DF_COLUMNS
        if column not in tracker_full_df.columns
    ]
    if missing_columns:
        raise ValueError(
            "tracker_full_df is missing columns needed for virtual packets: "
            f"{missing_columns}."
        )
    # Events of a packet must be contiguous and in order for replay
    if not tracker_full_df.index.is_monotonic_increasing:
        tracker_full_df = tracker_full_df.sort_index()

    packet_id, event_counts = np.unique(
        tracker_full_df.index.get_level_values("packet_id").to_numpy(),
        return_counts=True,
    )
    event_offsets = np.zeros(len(packet_id) + 1, dtype=np.int64)
    np.cumsum(event_counts, out=event_offsets[1:])
    return TrackerEventArrays(
        packet_id=packet_id.astype(np.int64),
        event_offsets=event_offsets,
        interaction_type=interaction_type_to_codes(
            tracker_full_df["interaction_type"]
        ),
        radius=tracker_full_df["radius"].to_numpy(dtype=np.float64),
        before_shell_id=tracker_full_df["before_shell_id"].to_numpy(
            dtype=np.int64
        ),
        after_shell_id=tracker_full_df["after_shell_id"].to_numpy(
            dtype=np.int64
        ),
        before_nu=tracker_full_df["before_nu"].to_numpy(dtype=np.float64),
        after_nu=tracker_full_df["after_nu"].to_numpy(dtype=np.float64),
        after_mu=tracker_full_df["after_mu"].to_numpy(dtype=np.float64),
        after_energy=tracker_full_df["after_energy"].to_numpy(dtype=np.float64),
        line_absorb_id=tracker_full_df["line_absorb_id"].to_numpy(
            dtype=np.int64
        ),
        line_emit_id=tracker_full_df["line_emit_id"].to_numpy(dtype=np.int64),
        next_line_id=tracker_full_df["next_line_id"].to_numpy(dtype=np.int64),
    )


@njit(**njit_dict_no_parallel)
def create_vpacket_collections(
    no_of_packets: int,
    spectrum_frequency_grid: npt.NDArray[np.float64],
    vpacket_spawn_start_frequency: float,
    vpacket_spawn_end_frequency: float,
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
    vpacket_spawn_start_frequency : float
        Lowest real-packet frequency that spawns virtual packets [Hz].
    vpacket_spawn_end_frequency : float
        Highest real-packet frequency that spawns virtual packets [Hz].
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
            vpacket_spawn_start_frequency,
            vpacket_spawn_end_frequency,
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
    tracker_events: TrackerEventArrays,
    vpacket_collections: TypedList[VPacketCollection],
    tau_russian: float,
    survival_probability: float,
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
    tracker_events : TrackerEventArrays
        Completed real-packet event histories, e.g. from
        ``tracker_full_df_to_event_arrays``.
    vpacket_collections : numba.typed.List
        Unused virtual packet collection for each real packet, in the order
        of ``tracker_events.packet_id``, e.g. from
        ``create_vpacket_collections``. Filled and finalized in place.
    tau_russian : float
        Optical depth above which virtual packets are Russian-rouletted.
    survival_probability : float
        Probability that a Russian-rouletted virtual packet survives.
    trace_vpacket_volley_function
        Mode-specific virtual-packet propagation function.

    Returns
    -------
    numpy.ndarray
        Virtual-packet energy histogram on ``spectrum_frequency_grid``.

    Raises
    ------
    ValueError
        If there is not exactly one unused collection per real packet.
    """
    no_of_packets = len(tracker_events.packet_id)
    # The collections are indexed without bounds checks and appended to in
    # place, so a mismatched or reused list (even one finalized while empty,
    # whose arrays are truncated below its length) would silently corrupt
    # the result
    if len(vpacket_collections) != no_of_packets:
        raise ValueError(
            "vpacket_collections must have one collection per real packet."
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
        packet_id = tracker_events.packet_id[packet_idx]
        vpacket_collection = vpacket_collections[packet_idx]
        # Offset the configured seed for deterministic parallel packet sampling.
        vpacket_seed = montecarlo_configuration.MONTECARLO_SEED + packet_id
        np.random.seed(vpacket_seed)

        for event_idx in range(
            tracker_events.event_offsets[packet_idx],
            tracker_events.event_offsets[packet_idx + 1],
        ):
            interaction_type = tracker_events.interaction_type[event_idx]
            is_initial_launch = (
                interaction_type == InteractionType.BOUNDARY
                and tracker_events.before_shell_id[event_idx] == -1
                and tracker_events.after_shell_id[event_idx] == 0
            )
            is_physical_interaction = (
                interaction_type == InteractionType.LINE
                or interaction_type == InteractionType.ESCATTERING
            )
            if not (is_initial_launch or is_physical_interaction):
                continue

            r_packet = RPacket(
                tracker_events.radius[event_idx],
                tracker_events.after_mu[event_idx],
                tracker_events.after_nu[event_idx],
                tracker_events.after_energy[event_idx],
                vpacket_seed,
                packet_id,
            )
            r_packet.current_shell_id = tracker_events.after_shell_id[event_idx]
            r_packet.next_line_id = tracker_events.next_line_id[event_idx]

            if is_initial_launch:
                last_interaction_in_nu = 0.0
                last_interaction_in_r = geometry_state_numba.r_inner[0]
                last_interaction_type = InteractionType.NO_INTERACTION
                last_interaction_in_id = -1
                last_interaction_out_id = -1
                last_interaction_shell_id = -1
            else:
                last_interaction_in_nu = tracker_events.before_nu[event_idx]
                last_interaction_in_r = tracker_events.radius[event_idx]
                last_interaction_type = interaction_type
                last_interaction_in_id = tracker_events.line_absorb_id[
                    event_idx
                ]
                last_interaction_out_id = tracker_events.line_emit_id[event_idx]
                last_interaction_shell_id = tracker_events.before_shell_id[
                    event_idx
                ]

            trace_vpacket_volley_function(
                r_packet,
                vpacket_collection,
                geometry_state_numba,
                time_explosion,
                opacity_state_numba,
                montecarlo_configuration.ENABLE_FULL_RELATIVITY,
                tau_russian,
                survival_probability,
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
