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

