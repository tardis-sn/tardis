import numba as nb
import numpy as np
from numba.experimental import jitclass


@jitclass
class OpacityStateNumbaIIP:
    """
    IIP-specific opacity state that extends the base OpacityStateNumba
    with absorbing Markov chain probabilities for faster macroatom interactions.
    """

    electron_density: nb.float64[:]  # type: ignore[misc]
    t_electrons: nb.float64[:]  # type: ignore[misc]
    line_list_nu: nb.float64[:]  # type: ignore[misc]
    tau_sobolev: nb.float64[:, :]  # type: ignore[misc]
    transition_probabilities: nb.float64[:, :]  # type: ignore[misc]
    line2macro_level_upper: nb.int64[:]  # type: ignore[misc]
    macro_block_references: nb.int64[:]  # type: ignore[misc]
    transition_type: nb.int64[:]  # type: ignore[misc]
    destination_level_id: nb.int64[:]  # type: ignore[misc]
    transition_line_id: nb.int64[:]  # type: ignore[misc]
    bf_threshold_list_nu: nb.float64[:]  # type: ignore[misc]
    p_fb_deactivation: nb.float64[:, :]  # type: ignore[misc]
    photo_ion_nu_threshold_mins: nb.float64[:]  # type: ignore[misc]
    photo_ion_nu_threshold_maxs: nb.float64[:]  # type: ignore[misc]
    photo_ion_block_references: nb.int64[:]  # type: ignore[misc]
    chi_bf: nb.float64[:, :]  # type: ignore[misc]
    x_sect: nb.float64[:]  # type: ignore[misc]
    phot_nus: nb.float64[:]  # type: ignore[misc]
    ff_opacity_factor: nb.float64[:]  # type: ignore[misc]
    emissivities: nb.float64[:, :]  # type: ignore[misc]
    photo_ion_activation_idx: nb.int64[:]  # type: ignore[misc]
    k_packet_idx: nb.int64  # type: ignore[misc]
    # IIP-specific field
    absorbing_markov_probabilities: nb.float64[:, :, :]  # type: ignore[misc]

    def __init__(
        self,
        electron_density: np.ndarray,
        t_electrons: np.ndarray,
        line_list_nu: np.ndarray,
        tau_sobolev: np.ndarray,
        transition_probabilities: np.ndarray,
        line2macro_level_upper: np.ndarray,
        macro_block_references: np.ndarray,
        transition_type: np.ndarray,
        destination_level_id: np.ndarray,
        transition_line_id: np.ndarray,
        bf_threshold_list_nu: np.ndarray,
        p_fb_deactivation: np.ndarray,
        photo_ion_nu_threshold_mins: np.ndarray,
        photo_ion_nu_threshold_maxs: np.ndarray,
        photo_ion_block_references: np.ndarray,
        chi_bf: np.ndarray,
        x_sect: np.ndarray,
        phot_nus: np.ndarray,
        ff_opacity_factor: np.ndarray,
        emissivities: np.ndarray,
        photo_ion_activation_idx: np.ndarray,
        k_packet_idx: int,
        absorbing_markov_probabilities: np.ndarray,
    ) -> None:
        """
        Initialize IIP-specific Numba-compatible opacity state for Monte Carlo transport.

        Parameters
        ----------
        electron_density : numpy.ndarray
            Electron density in each shell [cm^-3].
        t_electrons : numpy.ndarray
            Electron temperature in each shell [K].
        line_list_nu : numpy.ndarray
            Frequencies of spectral lines [Hz].
        tau_sobolev : numpy.ndarray
            Sobolev optical depths for line transitions.
        transition_probabilities : numpy.ndarray
            Probabilities for macro atom transitions.
        line2macro_level_upper : numpy.ndarray
            Mapping from line indices to macro atom upper levels.
        macro_block_references : numpy.ndarray
            Block references for macro atom data.
        transition_type : numpy.ndarray
            Type identifiers for transitions.
        destination_level_id : numpy.ndarray
            Destination level indices for transitions.
        transition_line_id : numpy.ndarray
            Line indices for transitions.
        bf_threshold_list_nu : numpy.ndarray
            Bound-free threshold frequencies [Hz].
        p_fb_deactivation : numpy.ndarray
            Free-bound deactivation probabilities.
        photo_ion_nu_threshold_mins : numpy.ndarray
            Minimum photoionization threshold frequencies [Hz].
        photo_ion_nu_threshold_maxs : numpy.ndarray
            Maximum photoionization threshold frequencies [Hz].
        photo_ion_block_references : numpy.ndarray
            Block references for photoionization data.
        chi_bf : numpy.ndarray
            Bound-free absorption coefficients.
        x_sect : numpy.ndarray
            Photoionization cross sections [cm^2].
        phot_nus : numpy.ndarray
            Photoionization frequencies [Hz].
        ff_opacity_factor : numpy.ndarray
            Free-free opacity factors.
        emissivities : numpy.ndarray
            Emission coefficients for bound-free transitions.
        photo_ion_activation_idx : numpy.ndarray
            Indices for photoionization activation.
        k_packet_idx : int
            Index for k-packet handling.
        absorbing_markov_probabilities : numpy.ndarray
            Matrix B: Absorbing probabilities of the Markov-chain macro atom.
            Shape: (n_shells, n_states, n_states). For each shell, contains the probability
            of being absorbed in each destination state when starting from each source state.
        """
        self.electron_density = electron_density
        self.t_electrons = t_electrons
        self.line_list_nu = line_list_nu
        self.tau_sobolev = tau_sobolev
        self.bf_threshold_list_nu = bf_threshold_list_nu

        #### Macro Atom transition probabilities
        self.transition_probabilities = transition_probabilities
        self.line2macro_level_upper = line2macro_level_upper

        self.macro_block_references = macro_block_references
        self.transition_type = transition_type

        # Destination level is not needed and/or generated for downbranch
        self.destination_level_id = destination_level_id
        self.transition_line_id = transition_line_id
        self.p_fb_deactivation = p_fb_deactivation

        # Continuum Opacity Data
        self.photo_ion_nu_threshold_mins = photo_ion_nu_threshold_mins
        self.photo_ion_nu_threshold_maxs = photo_ion_nu_threshold_maxs

        self.photo_ion_block_references = photo_ion_block_references
        self.chi_bf = chi_bf
        self.x_sect = x_sect
        self.phot_nus = phot_nus
        self.ff_opacity_factor = ff_opacity_factor
        self.emissivities = emissivities
        self.photo_ion_activation_idx = photo_ion_activation_idx
        self.k_packet_idx = k_packet_idx

        # IIP-specific field
        self.absorbing_markov_probabilities = absorbing_markov_probabilities

    def __getitem__(self, i: slice) -> "OpacityStateNumbaIIP":
        """Get a shell or slice of shells of the attributes of the opacity state.

        Parameters
        ----------
        i : slice
            Shell slice. Will fail if slice is int since class only supports array types.

        Returns
        -------
        OpacityStateNumbaIIP
            A shallow copy of the current instance with sliced data.
        """
        return OpacityStateNumbaIIP(
            self.electron_density[i],
            self.t_electrons[i],
            self.line_list_nu,
            self.tau_sobolev[:, i],
            self.transition_probabilities[:, i],
            self.line2macro_level_upper,
            self.macro_block_references,
            self.transition_type,
            self.destination_level_id,
            self.transition_line_id,
            self.bf_threshold_list_nu,
            self.p_fb_deactivation,
            self.photo_ion_nu_threshold_mins,
            self.photo_ion_nu_threshold_maxs,
            self.photo_ion_block_references,
            self.chi_bf,
            self.x_sect,
            self.phot_nus,
            self.ff_opacity_factor,
            self.emissivities,
            self.photo_ion_activation_idx,
            self.k_packet_idx,
            self.absorbing_markov_probabilities[:, i],
        )
