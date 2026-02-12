"""Continuum interaction estimators for Monte Carlo simulations.

This module contains estimators for continuum process tracking.
"""

import numba as nb
import numpy as np
from numba.experimental import jitclass
from numba.typed import List

from tardis.transport.montecarlo import njit_dict_no_parallel


@jitclass
class EstimatorsContinuum:
    """
    Estimators for continuum interaction processes.

    Attributes
    ----------
    photo_ion_estimator : numpy.ndarray
        Photoionization rate estimator.
    stim_recomb_estimator : numpy.ndarray
        Stimulated recombination rate estimator.
    bf_heating_estimator : numpy.ndarray
        Bound-free heating rate estimator.
    stim_recomb_cooling_estimator : numpy.ndarray
        Stimulated recombination cooling rate estimator.
    ff_heating_estimator : numpy.ndarray
        Free-free heating rate estimator.
    photo_ion_estimator_statistics : numpy.ndarray
        Photoionization statistics counter.
    """

    photo_ion_estimator: nb.float64[:, :]  # type: ignore[misc]
    stim_recomb_estimator: nb.float64[:, :]  # type: ignore[misc]
    bf_heating_estimator: nb.float64[:, :]  # type: ignore[misc]
    stim_recomb_cooling_estimator: nb.float64[:, :]  # type: ignore[misc]
    ff_heating_estimator: nb.float64[:]  # type: ignore[misc]
    photo_ion_estimator_statistics: nb.int64[:, :]  # type: ignore[misc]

    def __init__(
        self,
        photo_ion_estimator: np.ndarray,
        stim_recomb_estimator: np.ndarray,
        bf_heating_estimator: np.ndarray,
        stim_recomb_cooling_estimator: np.ndarray,
        ff_heating_estimator: np.ndarray,
        photo_ion_estimator_statistics: np.ndarray,
    ) -> None:
        """
        Initialize EstimatorsContinuum with pre-allocated arrays.

        Parameters
        ----------
        photo_ion_estimator
            Array for photoionization rates.
        stim_recomb_estimator
            Array for stimulated recombination rates.
        bf_heating_estimator
            Array for bound-free heating rates.
        stim_recomb_cooling_estimator
            Array for stimulated recombination cooling rates.
        ff_heating_estimator
            Array for free-free heating rates.
        photo_ion_estimator_statistics
            Array for photoionization statistics.
        """
        self.photo_ion_estimator = photo_ion_estimator
        self.stim_recomb_estimator = stim_recomb_estimator
        self.bf_heating_estimator = bf_heating_estimator
        self.stim_recomb_cooling_estimator = stim_recomb_cooling_estimator
        self.ff_heating_estimator = ff_heating_estimator
        self.photo_ion_estimator_statistics = photo_ion_estimator_statistics

    def increment(self, other: "EstimatorsContinuum") -> None:
        """
        Increment estimators with values from another instance.

        Parameters
        ----------
        other
            Another EstimatorsContinuum instance to add.
        """
        self.photo_ion_estimator += other.photo_ion_estimator
        self.stim_recomb_estimator += other.stim_recomb_estimator
        self.bf_heating_estimator += other.bf_heating_estimator
        self.stim_recomb_cooling_estimator += (
            other.stim_recomb_cooling_estimator
        )
        self.ff_heating_estimator += other.ff_heating_estimator
        self.photo_ion_estimator_statistics += (
            other.photo_ion_estimator_statistics
        )


@nb.njit(**njit_dict_no_parallel)
def init_estimators_continuum(
    n_levels_bf_species_by_n_cells_tuple: tuple[int, int], n_cells: int
) -> EstimatorsContinuum:
    """
    Factory function to create and initialize EstimatorsContinuum.

    Parameters
    ----------
    n_levels_bf_species_by_n_cells_tuple
        Shape tuple for bound-free transitions (n_levels_bf_species, n_cells)
    n_cells
        Number of cells in the simulation

    Returns
    -------
    EstimatorsContinuum
        Initialized estimators with zero-filled arrays
    """
    photo_ion_estimator = np.zeros(
        n_levels_bf_species_by_n_cells_tuple, dtype=np.float64
    )
    stim_recomb_estimator = np.zeros(
        n_levels_bf_species_by_n_cells_tuple, dtype=np.float64
    )
    bf_heating_estimator = np.zeros(
        n_levels_bf_species_by_n_cells_tuple, dtype=np.float64
    )
    stim_recomb_cooling_estimator = np.zeros(
        n_levels_bf_species_by_n_cells_tuple, dtype=np.float64
    )
    ff_heating_estimator = np.zeros(n_cells, dtype=np.float64)
    photo_ion_estimator_statistics = np.zeros(
        n_levels_bf_species_by_n_cells_tuple, dtype=np.int64
    )

    return EstimatorsContinuum(
        photo_ion_estimator,
        stim_recomb_estimator,
        bf_heating_estimator,
        stim_recomb_cooling_estimator,
        ff_heating_estimator,
        photo_ion_estimator_statistics,
    )


@nb.njit(**njit_dict_no_parallel)
def create_estimators_continuum_list(
    n_levels_bf_species_by_n_cells_tuple: tuple[int, int],
    n_cells: int,
    number: int,
) -> List[EstimatorsContinuum]:
    """
    Factory function to create a list of EstimatorsContinuum instances.

    Parameters
    ----------
    n_levels_bf_species_by_n_cells_tuple
        Shape tuple for bound-free transitions (n_levels_bf_species, n_cells)
    n_cells
        Number of cells in the simulation
    number
        Number of estimator instances to create

    Returns
    -------
    numba.typed.List[EstimatorsContinuum]
        Typed list of EstimatorsContinuum instances
    """
    estimator_list = List()

    for _ in range(number):
        estimator_list.append(
            init_estimators_continuum(
                n_levels_bf_species_by_n_cells_tuple, n_cells
            )
        )

    return estimator_list
