import numpy as np
from numba import float64, int64, njit
from numba.experimental import jitclass

from tardis.transport.montecarlo import (
    njit_dict_no_parallel,
)

packet_collection_spec = [
    ("initial_radii", float64[:]),
    ("initial_nus", float64[:]),
    ("initial_mus", float64[:]),
    ("initial_energies", float64[:]),
    ("packet_seeds", int64[:]),
    ("time_of_simulation", float64),
    ("radiation_field_luminosity", float64),  #
    ("output_nus", float64[:]),
    ("output_energies", float64[:]),
]


@jitclass(packet_collection_spec)
class PacketCollection:
    def __init__(
        self,
        initial_radii,
        initial_nus,
        initial_mus,
        initial_energies,
        packet_seeds,
        radiation_field_luminosity,
    ):
        self.initial_radii = initial_radii
        self.initial_nus = initial_nus
        self.initial_mus = initial_mus
        self.initial_energies = initial_energies
        self.packet_seeds = packet_seeds
        self.radiation_field_luminosity = radiation_field_luminosity
        self.time_of_simulation = (
            1 / radiation_field_luminosity
        )  # 1 erg / luminosity
        self.output_nus = np.ones_like(initial_radii, dtype=np.float64) * -99.0
        self.output_energies = (
            np.ones_like(initial_radii, dtype=np.float64) * -99.0
        )


@njit(**njit_dict_no_parallel)
def initialize_last_interaction_tracker(no_of_packets):
    last_line_interaction_in_ids = -1 * np.ones(no_of_packets, dtype=np.int64)
    last_line_interaction_out_ids = -1 * np.ones(no_of_packets, dtype=np.int64)
    last_line_interaction_shell_ids = -1 * np.ones(
        no_of_packets, dtype=np.int64
    )
    last_interaction_types = -1 * np.ones(no_of_packets, dtype=np.int64)
    last_interaction_in_nus = np.zeros(no_of_packets, dtype=np.float64)

    return LastInteractionTracker(
        last_interaction_types,
        last_interaction_in_nus,
        last_line_interaction_in_ids,
        last_line_interaction_out_ids,
        last_line_interaction_shell_ids,
    )


last_interaction_tracker_spec = [
    ("types", int64[:]),
    ("in_nus", float64[:]),
    ("in_ids", int64[:]),
    ("out_ids", int64[:]),
    ("shell_ids", int64[:]),
]


@jitclass(last_interaction_tracker_spec)
class LastInteractionTracker:
    def __init__(
        self,
        types,
        in_nus,
        in_ids,
        out_ids,
        shell_ids,
    ):
        self.types = types
        self.in_nus = in_nus
        self.in_ids = in_ids
        self.out_ids = out_ids
        self.shell_ids = shell_ids

    def update_last_interaction(self, r_packet, i):
        self.types[i] = r_packet.last_interaction_type
        self.in_nus[i] = r_packet.last_interaction_in_nu
        self.in_ids[i] = r_packet.last_line_interaction_in_id
        self.out_ids[i] = r_packet.last_line_interaction_out_id
        self.shell_ids[i] = r_packet.last_line_interaction_shell_id


vpacket_collection_spec = [
    ("source_rpacket_index", int64),
    ("spectrum_frequency", float64[:]),
    ("v_packet_spawn_start_frequency", float64),
    ("v_packet_spawn_end_frequency", float64),
    ("nus", float64[:]),
    ("energies", float64[:]),
    ("initial_mus", float64[:]),
    ("initial_rs", float64[:]),
    ("idx", int64),
    ("number_of_vpackets", int64),
    ("length", int64),
    ("last_interaction_in_nu", float64[:]),
    ("last_interaction_type", int64[:]),
    ("last_interaction_in_id", int64[:]),
    ("last_interaction_out_id", int64[:]),
    ("last_interaction_shell_id", int64[:]),
]


@jitclass(vpacket_collection_spec)
class VPacketCollection:
    def __init__(
        self,
        source_rpacket_index,
        spectrum_frequency,
        v_packet_spawn_start_frequency,
        v_packet_spawn_end_frequency,
        number_of_vpackets,
        temporary_v_packet_bins,
    ):
        self.spectrum_frequency = spectrum_frequency
        self.v_packet_spawn_start_frequency = v_packet_spawn_start_frequency
        self.v_packet_spawn_end_frequency = v_packet_spawn_end_frequency
        self.nus = np.empty(temporary_v_packet_bins, dtype=np.float64)
        self.energies = np.empty(temporary_v_packet_bins, dtype=np.float64)
        self.initial_mus = np.empty(temporary_v_packet_bins, dtype=np.float64)
        self.initial_rs = np.empty(temporary_v_packet_bins, dtype=np.float64)
        self.number_of_vpackets = number_of_vpackets
        self.last_interaction_in_nu = np.zeros(
            temporary_v_packet_bins, dtype=np.float64
        )
        self.last_interaction_type = -1 * np.ones(
            temporary_v_packet_bins, dtype=np.int64
        )
        self.last_interaction_in_id = -1 * np.ones(
            temporary_v_packet_bins, dtype=np.int64
        )
        self.last_interaction_out_id = -1 * np.ones(
            temporary_v_packet_bins, dtype=np.int64
        )
        self.last_interaction_shell_id = -1 * np.ones(
            temporary_v_packet_bins, dtype=np.int64
        )
        self.idx = 0
        self.source_rpacket_index = source_rpacket_index
        self.length = temporary_v_packet_bins

    def add_packet(
        self,
        nu,
        energy,
        initial_mu,
        initial_r,
        last_interaction_in_nu,
        last_interaction_type,
        last_interaction_in_id,
        last_interaction_out_id,
        last_interaction_shell_id,
    ):
        """
        Add a packet to the vpacket collection and potentially resizing the arrays.

        Parameters
        ----------
        nu : float
            Frequency of the packet.
        energy : float
            Energy of the packet.
        initial_mu : float
            Initial mu of the packet.
        initial_r : float
            Initial r of the packet.
        last_interaction_in_nu : float
            Frequency of the last interaction of the packet.
        last_interaction_type : int
            Type of the last interaction of the packet.
        last_interaction_in_id : int
            ID of the last interaction in the packet.
        last_interaction_out_id : int
            ID of the last interaction out of the packet.
        last_interaction_shell_id : int
            ID of the last interaction shell of the packet.

        Returns
        -------
        None

        """
        if self.idx >= self.length:
            temp_length = self.length * 2 + self.number_of_vpackets
            temp_nus = np.empty(temp_length, dtype=np.float64)
            temp_energies = np.empty(temp_length, dtype=np.float64)
            temp_initial_mus = np.empty(temp_length, dtype=np.float64)
            temp_initial_rs = np.empty(temp_length, dtype=np.float64)
            temp_last_interaction_in_nu = np.empty(
                temp_length, dtype=np.float64
            )
            temp_last_interaction_type = np.empty(temp_length, dtype=np.int64)
            temp_last_interaction_in_id = np.empty(temp_length, dtype=np.int64)
            temp_last_interaction_out_id = np.empty(temp_length, dtype=np.int64)
            temp_last_interaction_shell_id = np.empty(
                temp_length, dtype=np.int64
            )

            temp_nus[: self.length] = self.nus
            temp_energies[: self.length] = self.energies
            temp_initial_mus[: self.length] = self.initial_mus
            temp_initial_rs[: self.length] = self.initial_rs
            temp_last_interaction_in_nu[
                : self.length
            ] = self.last_interaction_in_nu
            temp_last_interaction_type[
                : self.length
            ] = self.last_interaction_type
            temp_last_interaction_in_id[
                : self.length
            ] = self.last_interaction_in_id
            temp_last_interaction_out_id[
                : self.length
            ] = self.last_interaction_out_id
            temp_last_interaction_shell_id[
                : self.length
            ] = self.last_interaction_shell_id

            self.nus = temp_nus
            self.energies = temp_energies
            self.initial_mus = temp_initial_mus
            self.initial_rs = temp_initial_rs
            self.last_interaction_in_nu = temp_last_interaction_in_nu
            self.last_interaction_type = temp_last_interaction_type
            self.last_interaction_in_id = temp_last_interaction_in_id
            self.last_interaction_out_id = temp_last_interaction_out_id
            self.last_interaction_shell_id = temp_last_interaction_shell_id
            self.length = temp_length

        self.nus[self.idx] = nu
        self.energies[self.idx] = energy
        self.initial_mus[self.idx] = initial_mu
        self.initial_rs[self.idx] = initial_r
        self.last_interaction_in_nu[self.idx] = last_interaction_in_nu
        self.last_interaction_type[self.idx] = last_interaction_type
        self.last_interaction_in_id[self.idx] = last_interaction_in_id
        self.last_interaction_out_id[self.idx] = last_interaction_out_id
        self.last_interaction_shell_id[self.idx] = last_interaction_shell_id
        self.idx += 1

    def finalize_arrays(self):
        """
        Finalize the arrays by truncating them based on the current index.

        Returns
        -------
        None

        """
        self.nus = self.nus[: self.idx]
        self.energies = self.energies[: self.idx]
        self.initial_mus = self.initial_mus[: self.idx]
        self.initial_rs = self.initial_rs[: self.idx]
        self.last_interaction_in_nu = self.last_interaction_in_nu[: self.idx]
        self.last_interaction_type = self.last_interaction_type[: self.idx]
        self.last_interaction_in_id = self.last_interaction_in_id[: self.idx]
        self.last_interaction_out_id = self.last_interaction_out_id[: self.idx]
        self.last_interaction_shell_id = self.last_interaction_shell_id[
            : self.idx
        ]


@njit(**njit_dict_no_parallel)
def consolidate_vpacket_tracker(
    vpacket_collections, spectrum_frequency, start_frequency, end_frequency
):
    """
    Consolidate the vpacket trackers from multiple collections into a single vpacket tracker.

    Parameters
    ----------
    vpacket_collections : List[VPacketCollection]
        List of vpacket collections to consolidate.
    spectrum_frequency : ndarray
        Array of spectrum frequencies.

    Returns
    -------
    VPacketCollection
        Consolidated vpacket tracker.

    """
    vpacket_tracker_length = 0
    for vpacket_collection in vpacket_collections:
        vpacket_tracker_length += vpacket_collection.idx

    vpacket_tracker = VPacketCollection(
        -1,
        spectrum_frequency,
        start_frequency,
        end_frequency,
        -1,
        vpacket_tracker_length,
    )
    current_start_vpacket_tracker_idx = 0
    for vpacket_collection in vpacket_collections:
        current_end_vpacket_tracker_idx = (
            current_start_vpacket_tracker_idx + vpacket_collection.idx
        )
        vpacket_tracker.nus[
            current_start_vpacket_tracker_idx:current_end_vpacket_tracker_idx
        ] = vpacket_collection.nus
        vpacket_tracker.energies[
            current_start_vpacket_tracker_idx:current_end_vpacket_tracker_idx
        ] = vpacket_collection.energies
        vpacket_tracker.initial_mus[
            current_start_vpacket_tracker_idx:current_end_vpacket_tracker_idx
        ] = vpacket_collection.initial_mus
        vpacket_tracker.initial_rs[
            current_start_vpacket_tracker_idx:current_end_vpacket_tracker_idx
        ] = vpacket_collection.initial_rs
        vpacket_tracker.last_interaction_in_nu[
            current_start_vpacket_tracker_idx:current_end_vpacket_tracker_idx
        ] = vpacket_collection.last_interaction_in_nu

        vpacket_tracker.last_interaction_type[
            current_start_vpacket_tracker_idx:current_end_vpacket_tracker_idx
        ] = vpacket_collection.last_interaction_type

        vpacket_tracker.last_interaction_in_id[
            current_start_vpacket_tracker_idx:current_end_vpacket_tracker_idx
        ] = vpacket_collection.last_interaction_in_id

        vpacket_tracker.last_interaction_out_id[
            current_start_vpacket_tracker_idx:current_end_vpacket_tracker_idx
        ] = vpacket_collection.last_interaction_out_id

        vpacket_tracker.last_interaction_shell_id[
            current_start_vpacket_tracker_idx:current_end_vpacket_tracker_idx
        ] = vpacket_collection.last_interaction_shell_id

        current_start_vpacket_tracker_idx = current_end_vpacket_tracker_idx
    return vpacket_tracker
