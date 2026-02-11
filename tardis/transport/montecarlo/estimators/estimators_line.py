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
    mean_intensity_blue : numpy.ndarray
        Blue-shifted mean intensity estimator per line and cell.
    energy_deposition_line : numpy.ndarray
        Energy deposition estimator per line and cell.
    """

    mean_intensity_blue: nb.float64[:, :]  # type: ignore[misc]
    energy_deposition_line: nb.float64[:, :]  # type: ignore[misc]

    def __init__(
        self,
        mean_intensity_blue: np.ndarray,
        energy_deposition_line: np.ndarray,
    ) -> None:
        """
        Initialize EstimatorsLine with pre-allocated arrays.

        Parameters
        ----------
        mean_intensity_blue
            Array for blue-shifted mean intensity per line and cell.
        energy_deposition_line
            Array for energy deposition per line and cell.
        """
        self.mean_intensity_blue = mean_intensity_blue
        self.energy_deposition_line = energy_deposition_line

    def increment(self, other: "EstimatorsLine") -> None:
        """
        Increment estimators with values from another instance.

        Parameters
        ----------
        other
            Another EstimatorsLine instance to add.
        """
        self.mean_intensity_blue += other.mean_intensity_blue
        self.energy_deposition_line += other.energy_deposition_line


@nb.njit(**njit_dict_no_parallel)
def init_estimators_line(
    tau_sobolev_shape: tuple[int, int],
) -> EstimatorsLine:
    """
    Factory function to create and initialize EstimatorsLine.

    Parameters
    ----------
    tau_sobolev_shape
        Shape of tau_sobolev array (n_lines, n_cells).

    Returns
    -------
    EstimatorsLine
        Initialized estimators with zero-filled arrays.
    """
    mean_intensity_blue = np.zeros(tau_sobolev_shape, dtype=np.float64)
    energy_deposition_line = np.zeros(tau_sobolev_shape, dtype=np.float64)

    return EstimatorsLine(mean_intensity_blue, energy_deposition_line)


@nb.njit(**njit_dict_no_parallel)
def create_estimators_line_list(
    tau_sobolev_shape: tuple[int, int], number: int
) -> List[EstimatorsLine]:
    """
    Factory function to create a list of EstimatorsLine instances.

    Parameters
    ----------
    tau_sobolev_shape
        Shape of tau_sobolev array (n_lines, n_cells).
    number
        Number of estimator instances to create.

    Returns
    -------
    numba.typed.List[EstimatorsLine]
        Typed list of EstimatorsLine instances.
    """
    estimator_list = List()

    for _ in range(number):
        estimator_list.append(init_estimators_line(tau_sobolev_shape))

    return estimator_list
