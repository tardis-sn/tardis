from typing import Self

import numba as nb
import numpy as np
import numpy.typing as npt
from numba.experimental import jitclass


@jitclass
class ContinuumOpacityStateNumba:
    """Array-backed continuum data used by IIP transport."""

    bf_threshold_list_nu: nb.float64[:]
    p_fb_deactivation: nb.float64[:, :]
    photo_ion_nu_threshold_mins: nb.float64[:]
    photo_ion_nu_threshold_maxs: nb.float64[:]
    photo_ion_block_references: nb.int64[:]
    chi_bf: nb.float64[:, :]
    x_sect: nb.float64[:]
    phot_nus: nb.float64[:]
    ff_opacity_factor: nb.float64[:]
    emissivities: nb.float64[:, :]
    photo_ion_activation_idx: nb.int64[:]
    k_packet_idx: nb.int64
    absorbing_markov_probabilities: nb.float64[:, :, :]
    macro_level_to_reduced_idx: nb.int64[:]
    reduced_to_macro_level_idx: nb.int64[:]

    def __init__(
        self,
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
        macro_level_to_reduced_idx: npt.NDArray[np.int64],
        reduced_to_macro_level_idx: npt.NDArray[np.int64],
    ) -> None:
        self.bf_threshold_list_nu = bf_threshold_list_nu
        self.p_fb_deactivation = p_fb_deactivation
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
        self.macro_level_to_reduced_idx = macro_level_to_reduced_idx
        self.reduced_to_macro_level_idx = reduced_to_macro_level_idx

    def __getitem__(self, i: slice) -> Self:
        """Return a shell-sliced view of this continuum state.

        Per-shell arrays are sliced along their shell axis. Per-continuum and
        per-frequency data, such as thresholds and cross-sections, are shared.
        """
        return ContinuumOpacityStateNumba(
            self.bf_threshold_list_nu,
            self.p_fb_deactivation[:, i],
            self.photo_ion_nu_threshold_mins,
            self.photo_ion_nu_threshold_maxs,
            self.photo_ion_block_references,
            self.chi_bf[:, i],
            self.x_sect,
            self.phot_nus,
            self.ff_opacity_factor[i],
            self.emissivities[:, i],
            self.photo_ion_activation_idx,
            self.k_packet_idx,
            self.absorbing_markov_probabilities[i],
            self.macro_level_to_reduced_idx,
            self.reduced_to_macro_level_idx,
        )


def create_empty_continuum_state_numba(
    n_shells: int,
) -> ContinuumOpacityStateNumba:
    """Create a continuum state with no bound-free continua and no free-free opacity.

    Transport without continuum processes carries this state so that every
    opacity state has the same Numba type. There are no bound-free thresholds,
    so no bound-free interpolation grid exists, and the zero free-free factors
    give zero free-free opacity in every shell.

    Parameters
    ----------
    n_shells : int
        Number of shells in the model.

    Returns
    -------
    ContinuumOpacityStateNumba
        Continuum state with zero-length bound-free arrays and zero free-free
        opacity factors.
    """
    return ContinuumOpacityStateNumba(
        np.zeros(0),
        np.zeros((0, n_shells)),
        np.zeros(0),
        np.zeros(0),
        np.zeros(1, dtype=np.int64),
        np.zeros((0, n_shells)),
        np.zeros(0),
        np.zeros(0),
        np.zeros(n_shells),
        np.zeros((0, n_shells)),
        np.zeros(0, dtype=np.int64),
        -1,
        np.zeros((n_shells, 0, 0)),
        np.zeros(0, dtype=np.int64),
        np.zeros(0, dtype=np.int64),
    )
