import numba as nb
import numpy as np
from numba import njit
from numba.experimental import jitclass
from numba.typed import List

from tardis.transport.montecarlo import njit_dict_no_parallel


@jitclass
class EstimatorsRadField:
    j_estimator: nb.float64[:]  # type: ignore[misc]
    nu_bar_estimator: nb.float64[:]  # type: ignore[misc]
    j_blue_estimator: nb.float64[:, :]  # type: ignore[misc]
    Edotlu_estimator: nb.float64[:, :]  # type: ignore[misc]

    def __init__(
        self,
        j_estimator: np.ndarray,
        nu_bar_estimator: np.ndarray,
        j_blue_estimator: np.ndarray,
        Edotlu_estimator: np.ndarray,
    ):
        self.j_estimator = j_estimator
        self.nu_bar_estimator = nu_bar_estimator
        self.j_blue_estimator = j_blue_estimator
        self.Edotlu_estimator = Edotlu_estimator

    def increment(self, other):
        """
        Increments each estimator with the corresponding estimator from another instance.

        Parameters
        ----------
        other : EstimatorsRadField
            Another instance of EstimatorsRadField.
        """
        self.j_estimator += other.j_estimator
        self.nu_bar_estimator += other.nu_bar_estimator
        self.j_blue_estimator += other.j_blue_estimator
        self.Edotlu_estimator += other.Edotlu_estimator

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
                EstimatorsRadField(
                    np.copy(self.j_estimator),
                    np.copy(self.nu_bar_estimator),
                    np.copy(self.j_blue_estimator),
                    np.copy(self.Edotlu_estimator),
                )
            )
        return estimator_list


@jitclass
class EstimatorsContinuum:
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
    ):
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
        other : EstimatorsContinuum
            Another instance of EstimatorsContinuum.
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
                EstimatorsContinuum(
                    np.copy(self.photo_ion_estimator),
                    np.copy(self.stim_recomb_estimator),
                    np.copy(self.bf_heating_estimator),
                    np.copy(self.stim_recomb_cooling_estimator),
                    np.copy(self.ff_heating_estimator),
                    np.copy(self.photo_ion_estimator_statistics),
                )
            )
        return estimator_list


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
def initialize_estimator_statistics(tau_sobolev_shape, gamma_shape):
    """
    Initializes the estimators used in the Monte Carlo simulation.

    Parameters
    ----------
    tau_sobolev_shape : tuple
        Shape of the array with the Sobolev optical depth.
    gamma_shape : tuple
        Shape of the array with the photoionization rate coefficients.

    Returns
    -------
    tuple
        A tuple containing (EstimatorsRadField, EstimatorsContinuum)
    """
    j_estimator = np.zeros(tau_sobolev_shape[1], dtype=np.float64)
    nu_bar_estimator = np.zeros(tau_sobolev_shape[1], dtype=np.float64)
    j_blue_estimator = np.zeros(tau_sobolev_shape, dtype=np.float64)
    Edotlu_estimator = np.zeros(tau_sobolev_shape, dtype=np.float64)

    photo_ion_estimator = np.zeros(gamma_shape, dtype=np.float64)
    stim_recomb_estimator = np.zeros(gamma_shape, dtype=np.float64)
    stim_recomb_cooling_estimator = np.zeros(gamma_shape, dtype=np.float64)
    bf_heating_estimator = np.zeros(gamma_shape, dtype=np.float64)
    ff_heating_estimator = np.zeros(gamma_shape[1], dtype=np.float64)
    photo_ion_estimator_statistics = np.zeros(gamma_shape, dtype=np.int64)

    radfield_estimators = EstimatorsRadField(
        j_estimator,
        nu_bar_estimator,
        j_blue_estimator,
        Edotlu_estimator,
    )

    continuum_estimators = EstimatorsContinuum(
        photo_ion_estimator,
        stim_recomb_estimator,
        bf_heating_estimator,
        stim_recomb_cooling_estimator,
        ff_heating_estimator,
        photo_ion_estimator_statistics,
    )

    return radfield_estimators, continuum_estimators
