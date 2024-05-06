import numpy as np
from numba import float64, int64
from numba.experimental import jitclass
from numba.typed import List


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
    Estimators
        The initialized estimators.

    Examples
    --------
    >>> tau_sobolev_shape = (10, 20)
    >>> gamma_shape = (5, 5)
    >>> initialize_estimators(tau_sobolev_shape, gamma_shape)
    <Estimators object at 0x...>
    """
    j_estimator = np.zeros(tau_sobolev_shape[1], dtype=np.float64)
    nu_bar_estimator = np.zeros(tau_sobolev_shape[1], dtype=np.float64)
    j_blue_estimator = np.zeros(tau_sobolev_shape)
    Edotlu_estimator = np.zeros(tau_sobolev_shape)

    photo_ion_estimator = np.zeros(gamma_shape, dtype=np.float64)
    stim_recomb_estimator = np.zeros(gamma_shape, dtype=np.float64)
    stim_recomb_cooling_estimator = np.zeros(gamma_shape, dtype=np.float64)
    bf_heating_estimator = np.zeros(gamma_shape, dtype=np.float64)

    stim_recomb_cooling_estimator = np.zeros(gamma_shape, dtype=np.float64)

    photo_ion_estimator_statistics = np.zeros(gamma_shape, dtype=np.int64)
    return RadiationFieldMCEstimators(
        j_estimator,
        nu_bar_estimator,
        j_blue_estimator,
        Edotlu_estimator,
        photo_ion_estimator,
        stim_recomb_estimator,
        bf_heating_estimator,
        stim_recomb_cooling_estimator,
        photo_ion_estimator_statistics,
    )


base_estimators_spec = [
    ("j_estimator", float64[:]),
    ("nu_bar_estimator", float64[:]),
    ("j_blue_estimator", float64[:, :]),
    ("Edotlu_estimator", float64[:, :]),
]

continuum_estimators_spec = [
    ("photo_ion_estimator", float64[:, :]),
    ("stim_recomb_estimator", float64[:, :]),
    ("bf_heating_estimator", float64[:, :]),
    ("stim_recomb_cooling_estimator", float64[:, :]),
    ("photo_ion_estimator_statistics", int64[:, :]),
]


@jitclass(base_estimators_spec + continuum_estimators_spec)
class RadiationFieldMCEstimators:
    def __init__(
        self,
        j_estimator,
        nu_bar_estimator,
        j_blue_estimator,
        Edotlu_estimator,
        photo_ion_estimator,
        stim_recomb_estimator,
        bf_heating_estimator,
        stim_recomb_cooling_estimator,
        photo_ion_estimator_statistics,
    ):
        self.j_estimator = j_estimator
        self.nu_bar_estimator = nu_bar_estimator
        self.j_blue_estimator = j_blue_estimator
        self.Edotlu_estimator = Edotlu_estimator
        self.photo_ion_estimator = photo_ion_estimator
        self.stim_recomb_estimator = stim_recomb_estimator
        self.bf_heating_estimator = bf_heating_estimator
        self.stim_recomb_cooling_estimator = stim_recomb_cooling_estimator
        self.photo_ion_estimator_statistics = photo_ion_estimator_statistics

    def increment(self, other):
        """
        Increments each estimator with the corresponding estimator from another instance of the class.

        Parameters
        ----------
        other : RadiationFieldMCEstimators
            Another instance of the RadiationFieldMCEstimators class.

        Returns
        -------
        None
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
        self.photo_ion_estimator_statistics += (
            other.photo_ion_estimator_statistics
        )

    def create_estimator_list(self, number):
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
                    np.copy(self.photo_ion_estimator_statistics),
                )
            )
        return estimator_list
