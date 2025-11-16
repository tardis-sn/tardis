import numba as nb
from numba.experimental import jitclass

from tardis.transport.montecarlo.packets.radiative_packet import InteractionType


@jitclass
class TrackerLastInteraction:
    r: nb.float64  # type: ignore[misc]
    nu: nb.float64  # type: ignore[misc]
    mu: nb.float64  # type: ignore[misc]
    energy: nb.float64  # type: ignore[misc]
    shell_id: nb.int64  # type: ignore[misc]
    interaction_type: nb.int64  # type: ignore[misc]

    # Interaction tracking (before/after)
    before_nu: nb.float64  # type: ignore[misc]
    before_mu: nb.float64  # type: ignore[misc]
    before_energy: nb.float64  # type: ignore[misc]
    after_nu: nb.float64  # type: ignore[misc]
    after_mu: nb.float64  # type: ignore[misc]
    after_energy: nb.float64  # type: ignore[misc]

    # Line interaction IDs (for line interactions only)
    interaction_line_absorb_id: nb.int64  # type: ignore[misc]
    interaction_line_emit_id: nb.int64  # type: ignore[misc]

    # Interaction counter
    interactions_count: nb.int64  # type: ignore[misc]

    """
    Numba JITCLASS for storing the last interaction the RPacket undergoes,
    with unified tracking for all interaction types.

    Parameters
    ----------
        r : float
            Radius of the shell where the RPacket is present
        nu : float
            Frequency of the RPacket
        energy : float
            Energy possessed by the RPacket
        shell_id : int
            Current Shell No in which the last interaction happened
        interaction_type: int
            Type of interaction the rpacket undergoes
        interaction_before_* : float
            Properties before interaction (nu, mu, energy)
        interaction_after_* : float
            Properties after interaction (nu, mu, energy)
        interaction_line_absorb_id : int
            Line ID for absorbed line interactions (-1 for non-line interactions)
        interaction_line_emit_id : int
            Line ID for emitted line interactions (-1 for non-line interactions)
        interactions_count : int
            Count of interactions
    """

    def __init__(self) -> None:
        """
        Initialize properties with default values.
        Float values are initialized with NaN for better data quality tracking.
        """
        self.r = -1.0
        self.nu = float("nan")
        self.mu = float("nan")  # MISSING INITIALIZATION - this was causing JIT vs non-JIT differences!
        self.energy = float("nan")
        self.shell_id = -1
        self.interaction_type = -1

        # Initialize interaction tracking
        self.before_nu = float("nan")
        self.before_mu = float("nan")
        self.before_energy = float("nan")
        self.after_nu = float("nan")
        self.after_mu = float("nan")
        self.after_energy = float("nan")

        # Initialize line interaction IDs
        self.interaction_line_absorb_id = -1
        self.interaction_line_emit_id = -1

        # Initialize counter
        self.interactions_count = 0

    def track_line_interaction_before(self, r_packet) -> None:
        """
        Track packet state before line interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet before line interaction.
        """
        self.before_nu = r_packet.nu
        self.before_mu = r_packet.mu
        self.before_energy = r_packet.energy
        # Track the line ID that will be absorbed
        self.interaction_line_absorb_id = r_packet.next_line_id

    def track_line_interaction_after(self, r_packet) -> None:
        """
        Track packet state after line interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet after line interaction.
        """
        self.after_nu = r_packet.nu
        self.after_mu = r_packet.mu
        self.after_energy = r_packet.energy
        # Track the line ID that was emitted (next_line_id - 1 after emission)
        self.interaction_line_emit_id = r_packet.next_line_id - 1
        self.interactions_count += 1
        # Update general tracking
        self.r = r_packet.r
        self.nu = r_packet.nu
        self.energy = r_packet.energy
        self.shell_id = r_packet.current_shell_id
        self.interaction_type = InteractionType.LINE  # Set interaction type directly

    def track_escattering_interaction_before(self, r_packet) -> None:
        """
        Track packet state before electron scattering interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet before electron scattering.
        """
        self.before_mu = r_packet.mu
        self.before_nu = r_packet.nu
        self.before_energy = r_packet.energy
        # Reset line IDs for non-line interactions
        self.interaction_line_absorb_id = -1
        self.interaction_line_emit_id = -1

    def track_escattering_interaction_after(self, r_packet) -> None:
        """
        Track packet state after electron scattering interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet after electron scattering.
        """
        self.after_mu = r_packet.mu
        self.after_nu = r_packet.nu
        self.after_energy = r_packet.energy
        self.interactions_count += 1
        # Update general tracking
        self.r = r_packet.r
        self.nu = r_packet.nu
        self.energy = r_packet.energy
        self.shell_id = r_packet.current_shell_id
        self.interaction_type = InteractionType.ESCATTERING  # Set interaction type directly

    def track_continuum_interaction_before(self, r_packet) -> None:
        """
        Track packet state before continuum process interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet before continuum process.
        """
        self.before_nu = r_packet.nu
        self.before_energy = r_packet.energy
        self.before_mu = r_packet.mu
        # Reset line IDs for non-line interactions
        self.interaction_line_absorb_id = -1
        self.interaction_line_emit_id = -1

    def track_continuum_interaction_after(self, r_packet) -> None:
        """
        Track packet state after continuum process interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet after continuum process.
        """
        self.after_nu = r_packet.nu
        self.after_energy = r_packet.energy
        self.after_mu = r_packet.mu
        self.interactions_count += 1
        # Update general tracking
        self.r = r_packet.r
        self.nu = r_packet.nu
        self.energy = r_packet.energy
        self.shell_id = r_packet.current_shell_id
        self.interaction_type = InteractionType.CONTINUUM_PROCESS  # Set interaction type directly

    def track_boundary_event(
        self, r_packet, from_shell_id=-1, to_shell_id=-1
    ) -> None:
        """
        Track packet state during boundary event.
        This method provides API compatibility with TrackerFull.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet during boundary event.
        from_shell_id : int, optional
            Shell ID the packet is leaving (default: -1).
        to_shell_id : int, optional
            Shell ID the packet is entering (default: -1).
        """
        # For last interaction tracker, boundary events don't count as interactions
        # But we update the position tracking
        self.r = r_packet.r
        self.shell_id = r_packet.current_shell_id

    def get_interaction_summary(self) -> nb.int64:  # type: ignore[misc]
        """
        Get summary of interaction count.

        Returns
        -------
        int
            Total interaction count.
        """
        return self.interactions_count


    def finalize(self) -> None:
        """
        Finalize tracker to match the common tracker API.

        TrackerLastInteraction stores only scalar values, so there is nothing
        to trim or consolidate. This method exists for interface compatibility
        with TrackerFull and is a no-op.
        """
        pass



