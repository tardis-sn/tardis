import numba as nb
import numpy as np
import numpy.typing as npt
from numba.experimental import jitclass


@jitclass
class OpacityStateNumbaIIP:
    """Store array-backed opacity and IIP interaction data for transport.

    This is a separate Numba ``jitclass`` because IIP transport requires the
    absorbing Markov-chain probability matrix in addition to the fields used
    by :class:`OpacityStateNumba`.
    """

    electron_density: nb.float64[:]  # type: ignore[misc]
    t_electrons: nb.float64[:]  # type: ignore[misc]
    line_list_nu: nb.float64[:]  # type: ignore[misc]
    tau_sobolev: nb.float64[:, :]  # type: ignore[misc]
    transition_probabilities: nb.float64[:, :]  # type: ignore[misc]
    line2macro_level_upper: nb.int64[:]  # type: ignore[misc]
    macro_block_edge_index: nb.int64[:]  # type: ignore[misc]
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
    absorbing_markov_probabilities: nb.float64[:, :, :]  # type: ignore[misc]

    def __init__(
        self,
        electron_density: npt.NDArray[np.float64],
        t_electrons: npt.NDArray[np.float64],
        line_list_nu: npt.NDArray[np.float64],
        tau_sobolev: npt.NDArray[np.float64],
        transition_probabilities: npt.NDArray[np.float64],
        line2macro_level_upper: npt.NDArray[np.int64],
        macro_block_edge_index: npt.NDArray[np.int64],
        transition_type: npt.NDArray[np.int64],
        destination_level_id: npt.NDArray[np.int64],
        transition_line_id: npt.NDArray[np.int64],
        bf_threshold_list_nu: npt.NDArray[np.float64],
        p_fb_deactivation: npt.NDArray[np.float64],
        photo_ion_nu_threshold_mins: npt.NDArray[np.float64],
        photo_ion_nu_threshold_maxs: npt.NDArray[np.float64],
        photo_ion_block_references: npt.NDArray[np.int64],
        chi_bf: npt.NDArray[np.float64],
        x_sect: npt.NDArray[np.float64],
        phot_nus: npt.NDArray[np.float64],
        ff_opacity_factor: npt.NDArray[np.float64],
        emissivities: npt.NDArray[np.float64],
        photo_ion_activation_idx: npt.NDArray[np.int64],
        k_packet_idx: int,
        absorbing_markov_probabilities: npt.NDArray[np.float64],
    ) -> None:
        """
        Initialize the Numba-compatible opacity and IIP interaction state.

        Parameters
        ----------
        electron_density
            Electron density in each shell [cm^-3].
        t_electrons
            Electron temperature in each shell [K].
        line_list_nu
            Frequencies of spectral lines [Hz].
        tau_sobolev
            Sobolev optical depths for line transitions.
        transition_probabilities
            Probabilities for macro atom transitions.
        line2macro_level_upper
            Mapping from line indices to macro atom upper levels.
        macro_block_edge_index
            Block references for macro atom data.
        transition_type
            Type identifiers for transitions.
        destination_level_id
            Destination level indices for transitions.
        transition_line_id
            Line indices for transitions.
        bf_threshold_list_nu
            Bound-free threshold frequencies [Hz].
        p_fb_deactivation
            Free-bound deactivation probabilities.
        photo_ion_nu_threshold_mins
            Minimum photoionization threshold frequencies [Hz].
        photo_ion_nu_threshold_maxs
            Maximum photoionization threshold frequencies [Hz].
        photo_ion_block_references
            Block references for photoionization data.
        chi_bf
            Bound-free absorption coefficients.
        x_sect
            Photoionization cross sections [cm^2].
        phot_nus
            Photoionization frequencies [Hz].
        ff_opacity_factor
            Free-free opacity factors.
        emissivities
            Emission coefficients for bound-free transitions.
        photo_ion_activation_idx
            Indices for photoionization activation.
        k_packet_idx
            Index for k-packet handling.
        absorbing_markov_probabilities
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

        self.macro_block_edge_index = macro_block_edge_index
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
        self.absorbing_markov_probabilities = absorbing_markov_probabilities

    def __getitem__(self, i: slice):
        """Return a shallow IIP state view for a contiguous shell slice.

        Parameters
        ----------
        i : slice
            Shell slice. Integer indexing is unsupported because the transport
            kernels require array-valued shell dimensions.

        Returns
        -------
        OpacityStateNumbaIIP
            A shallow state copy with the shell-dependent line, macro-atom,
            and absorbing-probability arrays sliced. Continuum arrays are not
            sliced and must remain compatible with the selected shell range.
        """
        return OpacityStateNumbaIIP(
            self.electron_density[i],
            self.t_electrons[i],
            self.line_list_nu,
            self.tau_sobolev[:, i],
            self.transition_probabilities[:, i],
            self.line2macro_level_upper,
            self.macro_block_edge_index,
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
            self.absorbing_markov_probabilities[i, :, :],
        )
