import numba as nb
import numpy as np
import pandas as pd
from numba import from_dtype, njit
from numba.experimental import jitclass
from pandas.api.types import CategoricalDtype

from tardis.transport.montecarlo.packets.radiative_packet import InteractionType

NO_INTERACTION_INT = int(InteractionType.NO_INTERACTION)

boundary_interaction_dtype = np.dtype(
    [
        ("event_id", "int64"),
        ("current_shell_id", "int64"),
        ("next_shell_id", "int64"),
    ]
)


@jitclass
class RPacketTracker:
    seed: nb.int64  # type: ignore[misc]
    status: nb.int64[:]  # type: ignore[misc]
    r: nb.float64[:]  # type: ignore[misc]
    nu: nb.float64[:]  # type: ignore[misc]
    mu: nb.float64[:]  # type: ignore[misc]
    energy: nb.float64[:]  # type: ignore[misc]
    shell_id: nb.int64[:]  # type: ignore[misc]
    interaction_type: nb.int64[:]  # type: ignore[misc]

    # LINE interaction tracking arrays (before/after for each interaction)
    line_before_nu: nb.float64[:]  # type: ignore[misc]
    line_before_mu: nb.float64[:]  # type: ignore[misc]
    line_before_id: nb.int64[:]  # type: ignore[misc]
    line_after_nu: nb.float64[:]  # type: ignore[misc]
    line_after_mu: nb.float64[:]  # type: ignore[misc]
    line_after_id: nb.int64[:]  # type: ignore[misc]

    # ESCATTERING interaction tracking arrays (before/after for each interaction)
    escat_before_mu: nb.float64[:]  # type: ignore[misc]
    escat_after_mu: nb.float64[:]  # type: ignore[misc]

    # CONTINUUM_PROCESS interaction tracking arrays (before/after for each interaction)
    continuum_before_nu: nb.float64[:]  # type: ignore[misc]
    continuum_before_energy: nb.float64[:]  # type: ignore[misc]
    continuum_before_mu: nb.float64[:]  # type: ignore[misc]
    continuum_after_nu: nb.float64[:]  # type: ignore[misc]
    continuum_after_energy: nb.float64[:]  # type: ignore[misc]
    continuum_after_mu: nb.float64[:]  # type: ignore[misc]

    boundary_interaction: from_dtype(boundary_interaction_dtype)[:]  # type: ignore[misc]
    num_interactions: nb.int64  # type: ignore[misc]
    boundary_interactions_index: nb.int64  # type: ignore[misc]
    event_id: nb.int64  # type: ignore[misc]
    extend_factor: nb.int64  # type: ignore[misc]
    """
    Numba JITCLASS for storing the information for each interaction a RPacket instance undergoes.

    Parameters
    ----------
        length : int
            Length of the initial array that is instantiated
        seed : int
            Seed for each RPacket
        index : int
            Index position of each RPacket
        status : int
            Current status of the RPacket as per interactions
        r : float
            Radius of the shell where the RPacket is present
        nu : float
            Frequency of the RPacket
        mu : float
            Cosine of the angle made by the direction of movement of the RPacket from its original direction
        energy : float
            Energy possessed by the RPacket at a particular shell
        shell_id : int
            Current Shell No in which the RPacket is present
        interaction_type: int
            Type of interaction the rpacket undergoes
        num_interactions : int
            Internal counter for the interactions that a particular RPacket undergoes
        extend_factor : int
            The factor by which to extend the properties array when the size limit is reached
    """

    def __init__(self, length: int) -> None:
        """
        Initialize the variables with default value
        """
        self.seed = np.int64(0)
        self.status = np.empty(length, dtype=np.int64)
        self.r = np.empty(length, dtype=np.float64)
        self.nu = np.empty(length, dtype=np.float64)
        self.mu = np.empty(length, dtype=np.float64)
        self.energy = np.empty(length, dtype=np.float64)
        self.shell_id = np.empty(length, dtype=np.int64)
        self.interaction_type = np.full(length, NO_INTERACTION_INT, dtype=np.int64)

        # Initialize LINE interaction tracking arrays
        self.line_before_nu = np.full(length, np.nan, dtype=np.float64)
        self.line_before_mu = np.full(length, np.nan, dtype=np.float64)
        self.line_before_id = np.full(length, -1, dtype=np.int64)
        self.line_after_nu = np.full(length, np.nan, dtype=np.float64)
        self.line_after_mu = np.full(length, np.nan, dtype=np.float64)
        self.line_after_id = np.full(length, -1, dtype=np.int64)

        # Initialize ESCATTERING interaction tracking arrays
        self.escat_before_mu = np.full(length, np.nan, dtype=np.float64)
        self.escat_after_mu = np.full(length, np.nan, dtype=np.float64)

        # Initialize CONTINUUM_PROCESS interaction tracking arrays
        self.continuum_before_nu = np.full(length, np.nan, dtype=np.float64)
        self.continuum_before_energy = np.full(length, np.nan, dtype=np.float64)
        self.continuum_before_mu = np.full(length, np.nan, dtype=np.float64)
        self.continuum_after_nu = np.full(length, np.nan, dtype=np.float64)
        self.continuum_after_energy = np.full(length, np.nan, dtype=np.float64)
        self.continuum_after_mu = np.full(length, np.nan, dtype=np.float64)

        self.boundary_interaction = np.empty(
            length,
            dtype=boundary_interaction_dtype,
        )
        self.num_interactions = 0
        self.boundary_interactions_index = 0
        self.event_id = 1
        self.extend_factor = 2

    def extend_array(self, array, array_length):
        temp_array = np.empty(
            array_length * self.extend_factor, dtype=array.dtype
        )
        temp_array[:array_length] = array
        return temp_array

    def extend_interaction_type_array(self, array, array_length):
        temp_array = np.full(
            array_length * self.extend_factor, NO_INTERACTION_INT, dtype=array.dtype
        )
        temp_array[:array_length] = array
        return temp_array

    def extend_float_array(self, array, array_length):
        temp_array = np.full(
            array_length * self.extend_factor, np.nan, dtype=array.dtype
        )
        temp_array[:array_length] = array
        return temp_array

    def extend_int_array(self, array, array_length):
        temp_array = np.full(
            array_length * self.extend_factor, -1, dtype=array.dtype
        )
        temp_array[:array_length] = array
        return temp_array

    def _extend_arrays_if_needed(self):
        """
        Extend all tracking arrays if needed based on current capacity.
        """
        if self.num_interactions >= self.status.size:
            self.status = self.extend_array(self.status, self.status.size)
            self.r = self.extend_array(self.r, self.r.size)
            self.nu = self.extend_array(self.nu, self.nu.size)
            self.mu = self.extend_array(self.mu, self.mu.size)
            self.energy = self.extend_array(self.energy, self.energy.size)
            self.shell_id = self.extend_array(self.shell_id, self.shell_id.size)
            self.interaction_type = self.extend_interaction_type_array(
                self.interaction_type, self.interaction_type.size
            )

            # Extend LINE interaction tracking arrays
            self.line_before_nu = self.extend_float_array(
                self.line_before_nu, self.line_before_nu.size
            )
            self.line_before_mu = self.extend_float_array(
                self.line_before_mu, self.line_before_mu.size
            )
            self.line_before_id = self.extend_int_array(
                self.line_before_id, self.line_before_id.size
            )
            self.line_after_nu = self.extend_float_array(
                self.line_after_nu, self.line_after_nu.size
            )
            self.line_after_mu = self.extend_float_array(
                self.line_after_mu, self.line_after_mu.size
            )
            self.line_after_id = self.extend_int_array(
                self.line_after_id, self.line_after_id.size
            )

            # Extend ESCATTERING interaction tracking arrays
            self.escat_before_mu = self.extend_float_array(
                self.escat_before_mu, self.escat_before_mu.size
            )
            self.escat_after_mu = self.extend_float_array(
                self.escat_after_mu, self.escat_after_mu.size
            )

            # Extend CONTINUUM_PROCESS interaction tracking arrays
            self.continuum_before_nu = self.extend_float_array(
                self.continuum_before_nu, self.continuum_before_nu.size
            )
            self.continuum_before_energy = self.extend_float_array(
                self.continuum_before_energy, self.continuum_before_energy.size
            )
            self.continuum_before_mu = self.extend_float_array(
                self.continuum_before_mu, self.continuum_before_mu.size
            )
            self.continuum_after_nu = self.extend_float_array(
                self.continuum_after_nu, self.continuum_after_nu.size
            )
            self.continuum_after_energy = self.extend_float_array(
                self.continuum_after_energy, self.continuum_after_energy.size
            )
            self.continuum_after_mu = self.extend_float_array(
                self.continuum_after_mu, self.continuum_after_mu.size
            )

    def _track_general_packet_state(self, r_packet):
        """
        Track general packet state properties.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet to track.
        """
        self.seed = r_packet.seed
        self.status[self.num_interactions] = r_packet.status
        self.r[self.num_interactions] = r_packet.r
        self.nu[self.num_interactions] = r_packet.nu
        self.mu[self.num_interactions] = r_packet.mu
        self.energy[self.num_interactions] = r_packet.energy
        self.shell_id[self.num_interactions] = r_packet.current_shell_id
        self.interaction_type[
            self.num_interactions
        ] = r_packet.last_interaction_type
        self.num_interactions += 1

    def track_line_interaction_before(self, r_packet) -> None:
        """
        Track packet state before line interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet before line interaction.
        """
        self._extend_arrays_if_needed()
        self.line_before_nu[self.num_interactions] = r_packet.nu
        self.line_before_mu[self.num_interactions] = r_packet.mu
        self.line_before_id[self.num_interactions] = r_packet.last_line_interaction_in_id

    def track_line_interaction_after(self, r_packet) -> None:
        """
        Track packet state after line interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet after line interaction.
        """
        self.line_after_nu[self.num_interactions] = r_packet.nu
        self.line_after_mu[self.num_interactions] = r_packet.mu
        self.line_after_id[self.num_interactions] = r_packet.last_line_interaction_out_id
        self._track_general_packet_state(r_packet)

    def track_escattering_interaction_before(self, r_packet) -> None:
        """
        Track packet state before electron scattering interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet before electron scattering.
        """
        self._extend_arrays_if_needed()
        self.escat_before_mu[self.num_interactions] = r_packet.mu

    def track_escattering_interaction_after(self, r_packet) -> None:
        """
        Track packet state after electron scattering interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet after electron scattering.
        """
        self.escat_after_mu[self.num_interactions] = r_packet.mu
        self._track_general_packet_state(r_packet)

    def track_continuum_interaction_before(self, r_packet) -> None:
        """
        Track packet state before continuum process interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet before continuum process.
        """
        self._extend_arrays_if_needed()
        self.continuum_before_nu[self.num_interactions] = r_packet.nu
        self.continuum_before_energy[self.num_interactions] = r_packet.energy
        self.continuum_before_mu[self.num_interactions] = r_packet.mu

    def track_continuum_interaction_after(self, r_packet) -> None:
        """
        Track packet state after continuum process interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet after continuum process.
        """
        self.continuum_after_nu[self.num_interactions] = r_packet.nu
        self.continuum_after_energy[self.num_interactions] = r_packet.energy
        self.continuum_after_mu[self.num_interactions] = r_packet.mu
        self._track_general_packet_state(r_packet)

    def track_boundary_interaction(self, current_shell_id, next_shell_id) -> None:
        """
        Track boundary interaction properties
        """
        if self.boundary_interactions_index >= self.boundary_interaction.size:
            self.boundary_interaction = self.extend_array(
                self.boundary_interaction,
                self.boundary_interaction.size,
            )

        self.boundary_interaction[self.boundary_interactions_index][
            "event_id"
        ] = self.event_id
        self.event_id += 1

        self.boundary_interaction[self.boundary_interactions_index][
            "current_shell_id"
        ] = current_shell_id

        self.boundary_interaction[self.boundary_interactions_index][
            "next_shell_id"
        ] = next_shell_id

        self.boundary_interactions_index += 1

    def track_boundary_crossing(self, r_packet, from_shell_id, to_shell_id) -> None:
        """
        Track boundary crossing - combines tracking packet state and boundary interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet during boundary crossing.
        from_shell_id : int
            Shell ID the packet is leaving.
        to_shell_id : int
            Shell ID the packet is entering.
        """
        # Extend arrays if needed and track general packet state
        self._extend_arrays_if_needed()
        self._track_general_packet_state(r_packet)
        # Track the specific boundary interaction
        self.track_boundary_interaction(from_shell_id, to_shell_id)

    def finalize_array(self):
        """
        Change the size of the array from length ( or multiple of length ) to
        the actual number of interactions
        """
        self.status = self.status[: self.num_interactions]
        self.r = self.r[: self.num_interactions]
        self.nu = self.nu[: self.num_interactions]
        self.mu = self.mu[: self.num_interactions]
        self.energy = self.energy[: self.num_interactions]
        self.shell_id = self.shell_id[: self.num_interactions]
        self.interaction_type = self.interaction_type[: self.num_interactions]

        # Finalize LINE interaction tracking arrays
        self.line_before_nu = self.line_before_nu[: self.num_interactions]
        self.line_before_mu = self.line_before_mu[: self.num_interactions]
        self.line_before_id = self.line_before_id[: self.num_interactions]
        self.line_after_nu = self.line_after_nu[: self.num_interactions]
        self.line_after_mu = self.line_after_mu[: self.num_interactions]
        self.line_after_id = self.line_after_id[: self.num_interactions]

        # Finalize ESCATTERING interaction tracking arrays
        self.escat_before_mu = self.escat_before_mu[: self.num_interactions]
        self.escat_after_mu = self.escat_after_mu[: self.num_interactions]

        # Finalize CONTINUUM_PROCESS interaction tracking arrays
        self.continuum_before_nu = self.continuum_before_nu[: self.num_interactions]
        self.continuum_before_energy = self.continuum_before_energy[: self.num_interactions]
        self.continuum_before_mu = self.continuum_before_mu[: self.num_interactions]
        self.continuum_after_nu = self.continuum_after_nu[: self.num_interactions]
        self.continuum_after_energy = self.continuum_after_energy[: self.num_interactions]
        self.continuum_after_mu = self.continuum_after_mu[: self.num_interactions]

        self.boundary_interaction = self.boundary_interaction[
            : self.boundary_interactions_index
        ]

    def track_line_interaction(self, interaction_idx, before_nu, before_mu, before_id, after_nu, after_mu, after_id):
        """Track line interaction details for a specific interaction."""
        self.line_before_nu[interaction_idx] = before_nu
        self.line_before_mu[interaction_idx] = before_mu
        self.line_before_id[interaction_idx] = before_id
        self.line_after_nu[interaction_idx] = after_nu
        self.line_after_mu[interaction_idx] = after_mu
        self.line_after_id[interaction_idx] = after_id

    def track_escattering_interaction(self, interaction_idx, before_mu, after_mu):
        """Track electron scattering interaction details for a specific interaction."""
        self.escat_before_mu[interaction_idx] = before_mu
        self.escat_after_mu[interaction_idx] = after_mu

    def track_continuum_interaction(self, interaction_idx, before_nu, before_energy, before_mu, after_nu, after_energy, after_mu):
        """Track continuum process interaction details for a specific interaction."""
        self.continuum_before_nu[interaction_idx] = before_nu
        self.continuum_before_energy[interaction_idx] = before_energy
        self.continuum_before_mu[interaction_idx] = before_mu
        self.continuum_after_nu[interaction_idx] = after_nu
        self.continuum_after_energy[interaction_idx] = after_energy
        self.continuum_after_mu[interaction_idx] = after_mu


def rpacket_trackers_to_dataframe(rpacket_trackers):
    """
    Convert a list of RPacketTracker instances to a DataFrame.

    This function extracts interaction data from each tracker and creates
    a pandas DataFrame with detailed interaction tracking similar to the
    LastInteractionTracker but for every interaction.

    Parameters
    ----------
    rpacket_trackers : numba.typed.typedlist.List
        List of RPacketTracker instances

    Returns
    -------
    pd.DataFrame
        DataFrame containing interaction data with columns for basic properties
        and interaction-specific details, indexed by packet_id and step
    """
    # Create categorical dtype from enum for better performance and type safety
    interaction_type_dtype = CategoricalDtype(
        categories=[member.name for member in InteractionType],
        ordered=False
    )

    # Collect all data arrays
    all_status = []
    all_seed = []
    all_r = []
    all_nu = []
    all_mu = []
    all_energy = []
    all_shell_id = []
    all_interaction_type = []
    all_packet_index = []
    all_step_index = []

    # Collect interaction-specific data
    all_line_before_nu = []
    all_line_before_mu = []
    all_line_before_id = []
    all_line_after_nu = []
    all_line_after_mu = []
    all_line_after_id = []
    all_escat_before_mu = []
    all_escat_after_mu = []
    all_continuum_before_nu = []
    all_continuum_before_energy = []
    all_continuum_before_mu = []
    all_continuum_after_nu = []
    all_continuum_after_energy = []
    all_continuum_after_mu = []

    for i, rpacket_tracker in enumerate(rpacket_trackers):
        n_interactions = len(rpacket_tracker.r)

        all_status.append(rpacket_tracker.status)
        all_seed.append(np.full(n_interactions, rpacket_tracker.seed))
        all_r.append(rpacket_tracker.r)
        all_nu.append(rpacket_tracker.nu)
        all_mu.append(rpacket_tracker.mu)
        all_energy.append(rpacket_tracker.energy)
        all_shell_id.append(rpacket_tracker.shell_id)
        all_interaction_type.append(rpacket_tracker.interaction_type)
        all_packet_index.append(np.full(n_interactions, i))  # Use i instead of rpacket_tracker.index
        all_step_index.append(np.arange(n_interactions))

        # Add interaction-specific data
        all_line_before_nu.append(rpacket_tracker.line_before_nu)
        all_line_before_mu.append(rpacket_tracker.line_before_mu)
        all_line_before_id.append(rpacket_tracker.line_before_id)
        all_line_after_nu.append(rpacket_tracker.line_after_nu)
        all_line_after_mu.append(rpacket_tracker.line_after_mu)
        all_line_after_id.append(rpacket_tracker.line_after_id)
        all_escat_before_mu.append(rpacket_tracker.escat_before_mu)
        all_escat_after_mu.append(rpacket_tracker.escat_after_mu)
        all_continuum_before_nu.append(rpacket_tracker.continuum_before_nu)
        all_continuum_before_energy.append(rpacket_tracker.continuum_before_energy)
        all_continuum_before_mu.append(rpacket_tracker.continuum_before_mu)
        all_continuum_after_nu.append(rpacket_tracker.continuum_after_nu)
        all_continuum_after_energy.append(rpacket_tracker.continuum_after_energy)
        all_continuum_after_mu.append(rpacket_tracker.continuum_after_mu)

    # Concatenate all arrays
    combined_interaction_type = np.concatenate(all_interaction_type)

    # Convert enum values to their string names and create categorical
    interaction_type_labels = [InteractionType(int_type).name for int_type in combined_interaction_type]
    interaction_type_categorical = pd.Categorical(
        interaction_type_labels,
        dtype=interaction_type_dtype
    )

    data = {
        "status": np.concatenate(all_status),
        "seed": np.concatenate(all_seed),
        "r": np.concatenate(all_r),
        "nu": np.concatenate(all_nu),
        "mu": np.concatenate(all_mu),
        "energy": np.concatenate(all_energy),
        "shell_id": np.concatenate(all_shell_id),
        "interaction_type": interaction_type_categorical,

        # LINE interaction columns
        "line_before_nu": np.concatenate(all_line_before_nu),
        "line_before_mu": np.concatenate(all_line_before_mu),
        "line_before_id": np.concatenate(all_line_before_id),
        "line_after_nu": np.concatenate(all_line_after_nu),
        "line_after_mu": np.concatenate(all_line_after_mu),
        "line_after_id": np.concatenate(all_line_after_id),

        # ESCATTERING interaction columns
        "escat_before_mu": np.concatenate(all_escat_before_mu),
        "escat_after_mu": np.concatenate(all_escat_after_mu),

        # CONTINUUM_PROCESS interaction columns
        "continuum_before_nu": np.concatenate(all_continuum_before_nu),
        "continuum_before_energy": np.concatenate(all_continuum_before_energy),
        "continuum_before_mu": np.concatenate(all_continuum_before_mu),
        "continuum_after_nu": np.concatenate(all_continuum_after_nu),
        "continuum_after_energy": np.concatenate(all_continuum_after_energy),
        "continuum_after_mu": np.concatenate(all_continuum_after_mu),
    }

    # Create multi-index with packet_id and step
    index_arrays = [np.concatenate(all_packet_index), np.concatenate(all_step_index)]
    multi_index = pd.MultiIndex.from_arrays(index_arrays, names=["packet_id", "step"])

    return pd.DataFrame(data, index=multi_index)
    combined_interaction_type = np.concatenate(all_interaction_type)

    data = {
        "status": np.concatenate(all_status),
        "seed": np.concatenate(all_seed),
        "r": np.concatenate(all_r),
        "nu": np.concatenate(all_nu),
        "mu": np.concatenate(all_mu),
        "energy": np.concatenate(all_energy),
        "shell_id": np.concatenate(all_shell_id),
        "interaction_type": combined_interaction_type,
    }

    # Create multi-index
    index_arrays = [np.concatenate(all_packet_index), np.concatenate(all_step_index)]
    multi_index = pd.MultiIndex.from_arrays(index_arrays, names=["index", "step"])

    return pd.DataFrame(data, index=multi_index)


@njit
def generate_rpacket_tracker_list(no_of_packets, length):
    """
    Parameters
    ----------
    no_of_packets : The count of RPackets that are sent in the ejecta
    length : initial length of the tracking array

    Returns
    -------
    A list containing RPacketTracker for each RPacket
    """
    # Create trackers list
    rpacket_trackers = []
    for i in range(no_of_packets):
        rpacket_trackers.append(RPacketTracker(length))
    return rpacket_trackers
