from typing import Self

import numba as nb
import numpy as np
import numpy.typing as npt
from numba.experimental import jitclass


@jitclass
class OpacityStateNumba:
    """Array-backed line opacity and macro-atom data for transport."""

    electron_density: nb.float64[:]
    t_electrons: nb.float64[:]
    line_list_nu: nb.float64[:]
    tau_sobolev: nb.float64[:, :]
    transition_probabilities: nb.float64[:, :]
    line2macro_level_upper: nb.int64[:]
    macro_block_edge_index: nb.int64[:]
    transition_type: nb.int64[:]
    destination_level_id: nb.int64[:]
    transition_line_id: nb.int64[:]

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
    ) -> None:
        self.electron_density = electron_density
        self.t_electrons = t_electrons
        self.line_list_nu = line_list_nu
        self.tau_sobolev = tau_sobolev
        self.transition_probabilities = transition_probabilities
        self.line2macro_level_upper = line2macro_level_upper
        self.macro_block_edge_index = macro_block_edge_index
        self.transition_type = transition_type
        self.destination_level_id = destination_level_id
        self.transition_line_id = transition_line_id

    def __getitem__(self, i: slice) -> Self:
        """Return a shell-sliced view of this opacity state."""
        return OpacityStateNumba(
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
        )
