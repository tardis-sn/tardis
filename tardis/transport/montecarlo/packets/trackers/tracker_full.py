import numba as nb
import numpy as np
from numba.experimental import jitclass

from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
)
from tardis.transport.montecarlo.packets.trackers.array_utils import (
    extend_array,
    extend_float_array,
    extend_int_array,
    extend_interaction_type_array,
)

NO_INTERACTION_INT = int(InteractionType.NO_INTERACTION)


@jitclass
class TrackerFull:
    """
    Numba JITCLASS for storing interaction information for RPacket instances.

    Tracks all packet events with before/after shell IDs and interaction data.

    Parameters
    ----------
    length : int
        Initial length of tracking arrays.

    Attributes
    ----------
    r : nb.float64[:]
        Radius of the RPacket at each event.
    shell_id : nb.int64[:]
        Shell ID where the RPacket is located at each event.
    before_shell_id : nb.int64[:]
        Shell ID before each event (same as shell_id except for boundary events).
    after_shell_id : nb.int64[:]
        Shell ID after each event (same as shell_id except for boundary events).
    interaction_type : nb.int64[:]
        Type of interaction undergone by the RPacket at each event.
    status : nb.int64[:]
        Status of the RPacket at each event.
    before_nu : nb.float64[:]
        Frequency before each interaction.
    before_mu : nb.float64[:]
        Cosine of angle before each interaction.
    before_energy : nb.float64[:]
        Energy before each interaction.
    after_nu : nb.float64[:]
        Frequency after each interaction.
    after_mu : nb.float64[:]
        Cosine of angle after each interaction.
    after_energy : nb.float64[:]
        Energy after each interaction.
    line_absorb_id : nb.int64[:]
        Line ID for absorbed line interactions.
    line_emit_id : nb.int64[:]
        Line ID for emitted line interactions.
    event_id : nb.int64
        Current event counter.
    extend_factor : nb.int64
        Factor by which to extend arrays when capacity is reached.
    """

    # All arrays are now "core" arrays
    radius: nb.float64[:]  # type: ignore[misc]
    before_shell_id: nb.int64[:]  # type: ignore[misc]
    after_shell_id: nb.int64[:]  # type: ignore[misc]
    interaction_type: nb.int64[:]  # type: ignore[misc]
    status: nb.int64[:]  # type: ignore[misc]
    before_nu: nb.float64[:]  # type: ignore[misc]
    before_mu: nb.float64[:]  # type: ignore[misc]
    before_energy: nb.float64[:]  # type: ignore[misc]
    after_nu: nb.float64[:]  # type: ignore[misc]
    after_mu: nb.float64[:]  # type: ignore[misc]
    after_energy: nb.float64[:]  # type: ignore[misc]
    line_absorb_id: nb.int64[:]  # type: ignore[misc]
    line_emit_id: nb.int64[:]  # type: ignore[misc]

    event_id: nb.int64  # type: ignore[misc]
    extend_factor: nb.int64  # type: ignore[misc]

    def __init__(self, length: int, extend_factor: int) -> None:
        """
        Initialize the RPacketTracker with arrays for tracking packet events.

        Parameters
        ----------
        length : int
            Initial length of the tracking arrays.
        """
        # All arrays have the same size now
        self.radius = np.full(length, np.nan, dtype=np.float64)
        self.before_shell_id = np.full(length, -99, dtype=np.int64)
        self.after_shell_id = np.full(length, -99, dtype=np.int64)
        self.interaction_type = np.full(
            length, NO_INTERACTION_INT, dtype=np.int64
        )
        self.status = np.full(length, -1, dtype=np.int64)
        self.before_nu = np.full(length, np.nan, dtype=np.float64)
        self.before_mu = np.full(length, np.nan, dtype=np.float64)
        self.before_energy = np.full(length, np.nan, dtype=np.float64)
        self.after_nu = np.full(length, np.nan, dtype=np.float64)
        self.after_mu = np.full(length, np.nan, dtype=np.float64)
        self.after_energy = np.full(length, np.nan, dtype=np.float64)
        self.line_absorb_id = np.full(length, -1, dtype=np.int64)
        self.line_emit_id = np.full(length, -1, dtype=np.int64)

        self.event_id = 0
        self.extend_factor = extend_factor

    @property
    def array_length(self) -> int:
        """
        Current capacity of tracking arrays.

        Returns
        -------
        int
            Current array capacity.
        """
        return self.radius.size

    def _extend_arrays(self) -> None:
        """Extend all tracking arrays."""
        new_length = self.array_length * self.extend_factor

        self.radius = extend_array(self.radius, new_length)
        self.before_shell_id = extend_array(self.before_shell_id, new_length)
        self.after_shell_id = extend_array(self.after_shell_id, new_length)
        self.interaction_type = extend_interaction_type_array(
            self.interaction_type, new_length
        )
        self.status = extend_array(self.status, new_length)
        self.before_nu = extend_float_array(self.before_nu, new_length)
        self.before_mu = extend_float_array(self.before_mu, new_length)
        self.before_energy = extend_float_array(self.before_energy, new_length)
        self.after_nu = extend_float_array(self.after_nu, new_length)
        self.after_mu = extend_float_array(self.after_mu, new_length)
        self.after_energy = extend_float_array(self.after_energy, new_length)
        self.line_absorb_id = extend_int_array(self.line_absorb_id, new_length)
        self.line_emit_id = extend_int_array(self.line_emit_id, new_length)

    def track_line_interaction_before(self, r_packet) -> None:
        """
        Track packet state before line interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet before line interaction.
        """
        # Extend arrays if needed
        if self.event_id >= self.array_length:
            self._extend_arrays()

        # Track all packet state
        self.radius[self.event_id] = r_packet.r
        self.before_shell_id[self.event_id] = r_packet.current_shell_id
        self.after_shell_id[self.event_id] = r_packet.current_shell_id
        self.interaction_type[self.event_id] = InteractionType.LINE
        self.status[self.event_id] = r_packet.status
        self.before_nu[self.event_id] = r_packet.nu
        self.before_mu[self.event_id] = r_packet.mu
        self.before_energy[self.event_id] = r_packet.energy
        self.line_absorb_id[self.event_id] = r_packet.next_line_id

    def track_line_interaction_after(self, r_packet) -> None:
        """
        Track packet state after line interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet after line interaction.
        """
        # Track after state for the current event
        self.after_nu[self.event_id] = r_packet.nu
        self.after_mu[self.event_id] = r_packet.mu
        self.after_energy[self.event_id] = r_packet.energy
        self.line_emit_id[self.event_id] = r_packet.next_line_id - 1

        # Increment event counter
        self.event_id += 1

    def track_escattering_interaction_before(self, r_packet) -> None:
        """
        Track packet state before electron scattering interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet before electron scattering.
        """
        # Extend arrays if needed
        if self.event_id >= self.array_length:
            self._extend_arrays()

        # Track all packet state
        self.radius[self.event_id] = r_packet.r
        self.before_shell_id[self.event_id] = r_packet.current_shell_id
        self.after_shell_id[self.event_id] = r_packet.current_shell_id
        self.interaction_type[self.event_id] = InteractionType.ESCATTERING
        self.status[self.event_id] = r_packet.status
        self.before_nu[self.event_id] = r_packet.nu
        self.before_mu[self.event_id] = r_packet.mu
        self.before_energy[self.event_id] = r_packet.energy

    def track_escattering_interaction_after(self, r_packet) -> None:
        """
        Track packet state after electron scattering interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet after electron scattering.
        """
        # Track after state for the current event
        self.after_nu[self.event_id] = r_packet.nu
        self.after_mu[self.event_id] = r_packet.mu
        self.after_energy[self.event_id] = r_packet.energy

        # Increment event counter
        self.event_id += 1

    def track_continuum_interaction_before(self, r_packet) -> None:
        """
        Track packet state before continuum process interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet before continuum process.
        """
        # Extend arrays if needed
        if self.event_id >= self.array_length:
            self._extend_arrays()

        # Track all packet state
        self.radius[self.event_id] = r_packet.r
        self.before_shell_id[self.event_id] = r_packet.current_shell_id
        self.after_shell_id[self.event_id] = r_packet.current_shell_id
        self.interaction_type[self.event_id] = InteractionType.CONTINUUM_PROCESS
        self.status[self.event_id] = r_packet.status
        self.before_nu[self.event_id] = r_packet.nu
        self.before_mu[self.event_id] = r_packet.mu
        self.before_energy[self.event_id] = r_packet.energy

    def track_continuum_interaction_after(self, r_packet) -> None:
        """
        Track packet state after continuum process interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet after continuum process.
        """
        # Track after state for the current event
        self.after_nu[self.event_id] = r_packet.nu
        self.after_mu[self.event_id] = r_packet.mu
        self.after_energy[self.event_id] = r_packet.energy

        # Increment event counter
        self.event_id += 1

    def track_boundary_event(
        self, r_packet, from_shell_id: int, to_shell_id: int
    ) -> None:
        """
        Track boundary event when packet crosses shell boundary.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet crossing the boundary.
        from_shell_id : int
            Shell ID the packet is leaving.
        to_shell_id : int
            Shell ID the packet is entering.
        """
        # Extend arrays if needed
        if self.event_id >= self.array_length:
            self._extend_arrays()

        # Track boundary event state
        self.radius[self.event_id] = r_packet.r
        self.before_shell_id[self.event_id] = from_shell_id
        self.after_shell_id[self.event_id] = to_shell_id
        self.interaction_type[self.event_id] = InteractionType.BOUNDARY
        self.status[self.event_id] = r_packet.status

        # Track before/after data for boundary events
        self.before_nu[self.event_id] = r_packet.nu
        self.before_mu[self.event_id] = r_packet.mu
        self.before_energy[self.event_id] = r_packet.energy
        # For boundary events, after values are the same as before (no change in packet properties)
        self.after_nu[self.event_id] = r_packet.nu
        self.after_mu[self.event_id] = r_packet.mu
        self.after_energy[self.event_id] = r_packet.energy

        # No line IDs for boundary events
        self.line_absorb_id[self.event_id] = -1
        self.line_emit_id[self.event_id] = -1

        self.event_id += 1

    def finalize(self) -> None:
        """
        Change the size of the arrays from length ( or multiple of length ) to
        the actual number of events.
        """
        # Finalize all arrays to actual event count
        self.radius = self.radius[: self.event_id]
        self.before_shell_id = self.before_shell_id[: self.event_id]
        self.after_shell_id = self.after_shell_id[: self.event_id]
        self.interaction_type = self.interaction_type[: self.event_id]
        self.status = self.status[: self.event_id]
        self.before_nu = self.before_nu[: self.event_id]
        self.before_mu = self.before_mu[: self.event_id]
        self.before_energy = self.before_energy[: self.event_id]
        self.after_nu = self.after_nu[: self.event_id]
        self.after_mu = self.after_mu[: self.event_id]
        self.after_energy = self.after_energy[: self.event_id]
        self.line_absorb_id = self.line_absorb_id[: self.event_id]
        self.line_emit_id = self.line_emit_id[: self.event_id]


