"""Bulk radiation field estimators for Monte Carlo simulations.

This module contains estimators for cell-level bulk radiation field properties.
"""

import numba as nb
import numpy as np
from numba.experimental import jitclass
from numba.typed import List

from tardis.transport.montecarlo import njit_dict_no_parallel


@jitclass
class EstimatorsBulk:
    """
    Estimators for bulk radiation field properties at the cell level.

    Attributes
    ----------
    mean_intensity_total : numpy.ndarray
        Mean intensity estimator per cell.
    mean_frequency : numpy.ndarray
        Mean frequency estimator per cell.
    """

    mean_intensity_total: nb.float64[:]  # type: ignore[misc]
    mean_frequency: nb.float64[:]  # type: ignore[misc]

    def __init__(
        self, mean_intensity_total: np.ndarray, mean_frequency: np.ndarray
    ) -> None:
        """
        Initialize EstimatorsBulk with pre-allocated arrays.

        Parameters
        ----------
        mean_intensity_total
            Array for mean intensity per cell.
        mean_frequency
            Array for mean frequency per cell.
        """
        self.mean_intensity_total = mean_intensity_total
        self.mean_frequency = mean_frequency

    def increment(self, other: "EstimatorsBulk") -> None:
        """
        Increment estimators with values from another instance.

        Parameters
        ----------
        other
            Another EstimatorsBulk instance to add.
        """
        self.mean_intensity_total += other.mean_intensity_total
        self.mean_frequency += other.mean_frequency


@nb.njit(**njit_dict_no_parallel)
def init_estimators_bulk(n_cells: int) -> EstimatorsBulk:
    """
    Factory function to create and initialize EstimatorsBulk.

    Parameters
    ----------
    n_cells
        Number of cells in the simulation.

    Returns
    -------
    EstimatorsBulk
        Initialized estimators with zero-filled arrays.
    """
    mean_intensity_total = np.zeros(n_cells, dtype=np.float64)
    mean_frequency = np.zeros(n_cells, dtype=np.float64)

    return EstimatorsBulk(mean_intensity_total, mean_frequency)


@nb.njit(**njit_dict_no_parallel)
def create_estimators_bulk_list(
    n_cells: int, number: int
) -> List[EstimatorsBulk]:
    """
    Factory function to create a list of EstimatorsBulk instances.

    Parameters
    ----------
    n_cells
        Number of cells in the simulation.
    number
        Number of estimator instances to create.

    Returns
    -------
    numba.typed.List[EstimatorsBulk]
        Typed list of EstimatorsBulk instances.
    """
    estimator_list = List()

    for _ in range(number):
        estimator_list.append(init_estimators_bulk(n_cells))

    return estimator_list
