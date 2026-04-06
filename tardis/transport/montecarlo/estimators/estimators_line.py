"""Line radiation field estimators for Monte Carlo simulations.

This module contains estimators for line-specific radiation field properties.
"""

import numba as nb
import numpy as np
from numba.experimental import jitclass
from numba.typed import List

from tardis.transport.montecarlo import njit_dict_no_parallel


@jitclass
class EstimatorsLine:
    """
    Estimators for line-specific radiation field properties.

    Attributes
    ----------
    mean_intensity_blueward : numpy.ndarray
        Monte Carlo estimator for the mean intensity in a line's extreme blue wing (Lucy 1999).
    energy_deposition_line_rate : numpy.ndarray
        Monte Carlo estimator for the rate per unit volume at which energy is removed
        from the radiation field by line excitations (Lucy 1999).
    """

    mean_intensity_blueward: nb.float64[:, :]  # type: ignore[misc]
    energy_deposition_line_rate: nb.float64[:, :]  # type: ignore[misc]

    def __init__(
        self,
        mean_intensity_blueward: np.ndarray,
        energy_deposition_line_rate: np.ndarray,
    ) -> None:
        """
        Initialize EstimatorsLine with pre-allocated arrays.

        Parameters
        ----------
        mean_intensity_blueward
            Monte Carlo estimator for the mean intensity in a line's extreme blue wing
        energy_deposition_line_rate
            Monte Carlo estimator for the rate per unit volume at which energy is removed
            from the radiation field by line excitations
        """
        self.mean_intensity_blueward = mean_intensity_blueward
        self.energy_deposition_line_rate = energy_deposition_line_rate

    def increment(self, other: "EstimatorsLine") -> None:
        """
        Increment estimators with values from another instance.

        Parameters
        ----------
        other
            Another EstimatorsLine instance to add.
        """
        self.mean_intensity_blueward += other.mean_intensity_blueward
        self.energy_deposition_line_rate += other.energy_deposition_line_rate


@nb.njit(**njit_dict_no_parallel)
def init_estimators_line(
    n_lines_by_n_cells_tuple: tuple[int, int],
) -> EstimatorsLine:
    """
    Factory function to create and initialize EstimatorsLine.

    Parameters
    ----------
    n_lines_by_n_cells_tuple
        Shape of tau_sobolev array (n_lines, n_cells).

    Returns
    -------
    Initialized estimators with zero-filled arrays.
    """
    mean_intensity_blueward = np.zeros(n_lines_by_n_cells_tuple, dtype=np.float64)
    energy_deposition_line_rate = np.zeros(n_lines_by_n_cells_tuple, dtype=np.float64)

    return EstimatorsLine(mean_intensity_blueward, energy_deposition_line_rate)


@nb.njit(**njit_dict_no_parallel)
def create_estimators_line_list(
    n_lines_by_n_cells_tuple: tuple[int, int], number: int
) -> List[EstimatorsLine]:
    """
    Factory function to create a list of EstimatorsLine instances.

    Parameters
    ----------
    n_lines_by_n_cells_tuple
        Shape tuple (n_lines, n_cells).
    number
        Number of estimator instances to create.

    Returns
    -------
    Typed list of EstimatorsLine instances.
    """
    estimator_list = List()

    for _ in range(number):
        estimator_list.append(init_estimators_line(n_lines_by_n_cells_tuple))

    return estimator_list
