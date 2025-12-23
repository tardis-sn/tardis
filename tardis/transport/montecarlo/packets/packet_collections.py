import numba as nb
import numpy as np
from numba import njit
from numba.experimental import jitclass

from tardis.transport.montecarlo import njit_dict_no_parallel
from tardis.transport.montecarlo.packets.radiative_packet import InteractionType

# Pre-calculate integer values for Numba compatibility
NO_INTERACTION_INT = int(InteractionType.NO_INTERACTION)


@jitclass
class PacketCollection:
    initial_radii: nb.float64[:]  # type: ignore[misc]
    initial_nus: nb.float64[:]  # type: ignore[misc]
    initial_mus: nb.float64[:]  # type: ignore[misc]
    initial_energies: nb.float64[:]  # type: ignore[misc]
    packet_seeds: nb.int64[:]  # type: ignore[misc]
    time_of_simulation: nb.float64  # type: ignore[misc]
    radiation_field_luminosity: nb.float64  # type: ignore[misc]
    output_nus: nb.float64[:]  # type: ignore[misc]
    output_energies: nb.float64[:]  # type: ignore[misc]

    def __init__(
        self,
        initial_radii: np.ndarray,
        initial_nus: np.ndarray,
        initial_mus: np.ndarray,
        initial_energies: np.ndarray,
        packet_seeds: np.ndarray,
        radiation_field_luminosity: float,
    ) -> None:
        """
        Initialize Numba-compatible packet collection for Monte Carlo transport.

        Parameters
        ----------
        initial_radii : numpy.ndarray
            Initial radii of packets [cm].
        initial_nus : numpy.ndarray
            Initial frequencies of packets [Hz].
        initial_mus : numpy.ndarray
            Initial directional cosines of packets.
        initial_energies : numpy.ndarray
            Initial energies of packets [erg].
        packet_seeds : numpy.ndarray
            Random number seeds for packets.
        radiation_field_luminosity : float
            Luminosity of the radiation field [erg/s].
        """
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

    @property
    def number_of_packets(self) -> int:
        """
        Get the number of packets in the collection.

        Returns
        -------
        int
            Number of packets.
        """
        return len(self.initial_radii)


@njit(**njit_dict_no_parallel)
def initialize_last_interaction_tracker(no_of_packets):
    last_line_interaction_in_ids = -1 * np.ones(no_of_packets, dtype=np.int64)
    last_line_interaction_out_ids = -1 * np.ones(no_of_packets, dtype=np.int64)
    last_line_interaction_shell_ids = -1 * np.ones(
        no_of_packets, dtype=np.int64
    )
    last_interaction_types = NO_INTERACTION_INT * np.ones(
        no_of_packets, dtype=np.int64
    )
    last_interaction_in_nus = np.zeros(no_of_packets, dtype=np.float64)
    last_interaction_in_rs = np.zeros(no_of_packets, dtype=np.float64)

    return LastInteractionTracker(
        last_interaction_types,
        last_interaction_in_nus,
        last_interaction_in_rs,
        last_line_interaction_in_ids,
        last_line_interaction_out_ids,
        last_line_interaction_shell_ids,
    )



@jitclass
class VPacketCollection:
    source_rpacket_index: nb.int64  # type: ignore[misc]
    spectrum_frequency_grid: nb.float64[:]  # type: ignore[misc]
    v_packet_spawn_start_frequency: nb.float64  # type: ignore[misc]
    v_packet_spawn_end_frequency: nb.float64  # type: ignore[misc]
    nus: nb.float64[:]  # type: ignore[misc]
    energies: nb.float64[:]  # type: ignore[misc]
    initial_mus: nb.float64[:]  # type: ignore[misc]
    initial_rs: nb.float64[:]  # type: ignore[misc]
    idx: nb.int64  # type: ignore[misc]
    number_of_vpackets: nb.int64  # type: ignore[misc]
    length: nb.int64  # type: ignore[misc]
    last_interaction_in_nu: nb.float64[:]  # type: ignore[misc]
    last_interaction_in_r: nb.float64[:]  # type: ignore[misc]
    last_interaction_type: nb.int64[:]  # type: ignore[misc]
    last_interaction_in_id: nb.int64[:]  # type: ignore[misc]
    last_interaction_out_id: nb.int64[:]  # type: ignore[misc]
    last_interaction_shell_id: nb.int64[:]  # type: ignore[misc]

    def __init__(
        self,
        source_rpacket_index: int,
        spectrum_frequency_grid: np.ndarray,
        v_packet_spawn_start_frequency: float,
        v_packet_spawn_end_frequency: float,
        number_of_vpackets: int,
        temporary_v_packet_bins: int,
    ) -> None:
        """
        Initialize virtual packet collection for Monte Carlo transport.

        Parameters
        ----------
        source_rpacket_index : int
            Index of the source R-packet.
        spectrum_frequency_grid : numpy.ndarray
            Frequency grid for spectrum calculation [Hz].
        v_packet_spawn_start_frequency : float
            Start frequency for virtual packet spawning [Hz].
        v_packet_spawn_end_frequency : float
            End frequency for virtual packet spawning [Hz].
        number_of_vpackets : int
            Number of virtual packets to generate.
        temporary_v_packet_bins : int
            Initial size of temporary storage arrays.
        """
        self.spectrum_frequency_grid = spectrum_frequency_grid
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
        self.last_interaction_in_r = np.zeros(
            temporary_v_packet_bins, dtype=np.float64
        )
        self.last_interaction_type = NO_INTERACTION_INT * np.ones(
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
        nu: float,
        energy: float,
        initial_mu: float,
        initial_r: float,
        last_interaction_in_nu: float,
        last_interaction_in_r: float,
        last_interaction_type: int,
        last_interaction_in_id: int,
        last_interaction_out_id: int,
        last_interaction_shell_id: int,
    ) -> None:
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
        last_interaction_in_r : float
            Radius of the last interaction of the packet.
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
            temp_last_interaction_in_r = np.empty(temp_length, dtype=np.float64)
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
            temp_last_interaction_in_nu[: self.length] = (
                self.last_interaction_in_nu
            )
            temp_last_interaction_in_r[: self.length] = (
                self.last_interaction_in_r
            )
            temp_last_interaction_type[: self.length] = (
                self.last_interaction_type
            )
            temp_last_interaction_in_id[: self.length] = (
                self.last_interaction_in_id
            )
            temp_last_interaction_out_id[: self.length] = (
                self.last_interaction_out_id
            )
            temp_last_interaction_shell_id[: self.length] = (
                self.last_interaction_shell_id
            )

            self.nus = temp_nus
            self.energies = temp_energies
            self.initial_mus = temp_initial_mus
            self.initial_rs = temp_initial_rs
            self.last_interaction_in_nu = temp_last_interaction_in_nu
            self.last_interaction_in_r = temp_last_interaction_in_r
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
        self.last_interaction_in_r[self.idx] = last_interaction_in_r
        self.last_interaction_type[self.idx] = last_interaction_type
        self.last_interaction_in_id[self.idx] = last_interaction_in_id
        self.last_interaction_out_id[self.idx] = last_interaction_out_id
        self.last_interaction_shell_id[self.idx] = last_interaction_shell_id
        self.idx += 1

    def finalize_arrays(self) -> None:
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
        self.last_interaction_in_r = self.last_interaction_in_r[: self.idx]
        self.last_interaction_type = self.last_interaction_type[: self.idx]
        self.last_interaction_in_id = self.last_interaction_in_id[: self.idx]
        self.last_interaction_out_id = self.last_interaction_out_id[: self.idx]
        self.last_interaction_shell_id = self.last_interaction_shell_id[
            : self.idx
        ]


@njit(**njit_dict_no_parallel)
def consolidate_vpacket_tracker(
    vpacket_collections,
    spectrum_frequency_grid: np.ndarray,
    start_frequency: float,
    end_frequency: float,
) -> "VPacketCollection":
    """
    Consolidate the vpacket trackers from multiple collections into a single vpacket tracker.

    Parameters
    ----------
    vpacket_collections : List[VPacketCollection]
        List of vpacket collections to consolidate.
    spectrum_frequency_grid : ndarray
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
        spectrum_frequency_grid,
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
        vpacket_tracker.last_interaction_in_r[
            current_start_vpacket_tracker_idx:current_end_vpacket_tracker_idx
        ] = vpacket_collection.last_interaction_in_r

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
