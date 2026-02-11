import numba as nb
import numpy as np
from numba import njit
from numba.experimental import jitclass
from numba.typed import List

from tardis.transport.montecarlo import njit_dict_no_parallel
from tardis.transport.montecarlo.estimators.estimators_bulk import (
    EstimatorsBulk,
    init_estimators_bulk,
)
from tardis.transport.montecarlo.estimators.estimators_continuum import (
    EstimatorsContinuum,
    init_estimators_continuum,
)
from tardis.transport.montecarlo.estimators.estimators_line import (
    EstimatorsLine,
    init_estimators_line,
)


@jitclass
class RadiationFieldMCEstimators:
    j_estimator: nb.float64[:]  # type: ignore[misc]
    nu_bar_estimator: nb.float64[:]  # type: ignore[misc]
    j_blue_estimator: nb.float64[:, :]  # type: ignore[misc]
    Edotlu_estimator: nb.float64[:, :]  # type: ignore[misc]
    photo_ion_estimator: nb.float64[:, :]  # type: ignore[misc]
    stim_recomb_estimator: nb.float64[:, :]  # type: ignore[misc]
    bf_heating_estimator: nb.float64[:, :]  # type: ignore[misc]
    stim_recomb_cooling_estimator: nb.float64[:, :]  # type: ignore[misc]
    ff_heating_estimator: nb.float64[:]  # type: ignore[misc]
    photo_ion_estimator_statistics: nb.int64[:, :]  # type: ignore[misc]

    def __init__(
        self,
        j_estimator: np.ndarray,
        nu_bar_estimator: np.ndarray,
        j_blue_estimator: np.ndarray,
        Edotlu_estimator: np.ndarray,
        photo_ion_estimator: np.ndarray,
        stim_recomb_estimator: np.ndarray,
        bf_heating_estimator: np.ndarray,
        stim_recomb_cooling_estimator: np.ndarray,
        ff_heating_estimator: np.ndarray,
        photo_ion_estimator_statistics: np.ndarray,
    ):
        self.j_estimator = j_estimator
        self.nu_bar_estimator = nu_bar_estimator
        self.j_blue_estimator = j_blue_estimator
        self.Edotlu_estimator = Edotlu_estimator
        self.photo_ion_estimator = photo_ion_estimator
        self.stim_recomb_estimator = stim_recomb_estimator
        self.bf_heating_estimator = bf_heating_estimator
        self.stim_recomb_cooling_estimator = stim_recomb_cooling_estimator
        self.ff_heating_estimator = ff_heating_estimator
        self.photo_ion_estimator_statistics = photo_ion_estimator_statistics

    def increment(self, other):
        """
        Increments each estimator with the corresponding estimator from another instance.

        Parameters
        ----------
        other : RadiationFieldMCEstimators
            Another instance of RadiationFieldMCEstimators.
        """
        self.j_estimator += other.j_estimator
        self.nu_bar_estimator += other.nu_bar_estimator
        self.j_blue_estimator += other.j_blue_estimator
        self.Edotlu_estimator += other.Edotlu_estimator
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

    def create_estimator_list(self, number: int):
        """
        Create a list of estimator copies for parallel processing.

        Parameters
        ----------
        number : int
            Number of estimator copies to create.
        """
        estimator_list = List()

        for i in range(number):
            estimator_list.append(
                RadiationFieldMCEstimators(
                    np.copy(self.j_estimator),
                    np.copy(self.nu_bar_estimator),
                    np.copy(self.j_blue_estimator),
                    np.copy(self.Edotlu_estimator),
                    np.copy(self.photo_ion_estimator),
                    np.copy(self.stim_recomb_estimator),
                    np.copy(self.bf_heating_estimator),
                    np.copy(self.stim_recomb_cooling_estimator),
                    np.copy(self.ff_heating_estimator),
                    np.copy(self.photo_ion_estimator_statistics),
                )
            )
        return estimator_list


@njit(**njit_dict_no_parallel)
def initialize_estimator_statistics(
    tau_sobolev_shape: tuple[int, int],
    gamma_shape: tuple[int, int],
    n_cells: int,
) -> tuple[EstimatorsBulk, EstimatorsLine, EstimatorsContinuum]:
    """
    Initialize Monte Carlo estimators for bulk, line, and continuum radiation field.

    Creates separate estimator objects for cell-level bulk estimators,
    line-level radiation estimators, and continuum interaction estimators.

    Parameters
    ----------
    tau_sobolev_shape
        Shape of the Sobolev optical depth array (n_lines, n_cells).
    gamma_shape
        Shape of the photoionization rate coefficient array (n_levels, n_cells).
    n_cells
        Number of spatial cells.

    Returns
    -------
        Tuple of (EstimatorsBulk, EstimatorsLine, EstimatorsContinuum).
    """
    estimators_bulk = init_estimators_bulk(n_cells)
    estimators_line = init_estimators_line(tau_sobolev_shape)
    estimators_continuum = init_estimators_continuum(gamma_shape, n_cells)

    return estimators_bulk, estimators_line, estimators_continuum
