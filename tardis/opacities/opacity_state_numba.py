from typing import Self

import numba as nb
import numpy as np
import numpy.typing as npt
from numba.experimental import jitclass

from tardis.opacities.continuum.continuum_state_numba import (
    ContinuumOpacityStateNumba,
)


@jitclass
class OpacityStateNumba:
    """Array-backed line, macro-atom and continuum opacity data for transport.

    Attributes
    ----------
    sobolev_optical_depth_coefficient : numpy.ndarray
        Velocity-gradient-independent coefficient for each line and shell [s^-1].
    continuum : ContinuumOpacityStateNumba
        Continuum data; empty when continuum processes are disabled.
    """

    electron_density: nb.float64[:]
    t_electrons: nb.float64[:]
    line_list_nu: nb.float64[:]
    tau_sobolev: nb.float64[:, :]
    sobolev_optical_depth_coefficient: nb.float64[:, :]
    transition_probabilities: nb.float64[:, :]
    line2macro_level_upper: nb.int64[:]
    macro_block_edge_index: nb.int64[:]
    transition_type: nb.int64[:]
    destination_level_id: nb.int64[:]
    transition_line_id: nb.int64[:]
    continuum: ContinuumOpacityStateNumba.class_type.instance_type

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
        continuum: ContinuumOpacityStateNumba,
    ) -> None:
        self.electron_density = electron_density
        self.t_electrons = t_electrons
        self.line_list_nu = line_list_nu
        self.tau_sobolev = tau_sobolev
        self.sobolev_optical_depth_coefficient = np.zeros_like(tau_sobolev)
        self.transition_probabilities = transition_probabilities
        self.line2macro_level_upper = line2macro_level_upper
        self.macro_block_edge_index = macro_block_edge_index
        self.transition_type = transition_type
        self.destination_level_id = destination_level_id
        self.transition_line_id = transition_line_id
        self.continuum = continuum

    def __getitem__(self, i: slice) -> Self:
        """Return a shell-sliced view of this opacity state, including its continuum."""
        sliced_opacity_state = OpacityStateNumba(
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
            self.continuum[i],
        )
        sliced_opacity_state.sobolev_optical_depth_coefficient = (
            self.sobolev_optical_depth_coefficient[:, i]
        )
        return sliced_opacity_state
