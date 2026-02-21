from dataclasses import dataclass

import pandas as pd

import tardis.constants as const
from tardis.io.atom_data import AtomData
from tardis.transport.montecarlo.estimators.estimators_continuum import (
    EstimatorsContinuum,
)
from tardis.transport.montecarlo.estimators.util import (
    bound_free_estimator_array2frame,
)

H = const.h.cgs.value


class MCContinuumPropertiesSolver:
    def __init__(
        self,
        atom_data: AtomData,
    ):
        self.atom_data = atom_data

    def solve(
        self,
        estimators_continuum: EstimatorsContinuum,
        time_simulation: float,
        volume: float,
    ) -> "ContinuumProperties":
        """
        Solve for the continuum properties.

        Parameters
        ----------
        estimators_continuum
            The Monte Carlo estimators for the continuum radiation field.
        time_simulation
            The simulation time.
        volume
            The volume of the cells.

        Returns
        -------
        The calculated continuum properties.
        """
        photo_ion_norm_factor = (time_simulation * volume * H) ** -1

        photo_ionization_rate_coefficient = bound_free_estimator_array2frame(
            estimators_continuum.photo_ion_estimator,
            self.atom_data.level2continuum_edge_idx,
        )
        photo_ionization_rate_coefficient *= photo_ion_norm_factor

        stimulated_recomb_rate_factor = bound_free_estimator_array2frame(
            estimators_continuum.stim_recomb_estimator,
            self.atom_data.level2continuum_edge_idx,
        )
        stimulated_recomb_rate_factor *= photo_ion_norm_factor

        return ContinuumProperties(
            stimulated_recomb_rate_factor, photo_ionization_rate_coefficient
        )


@dataclass
class ContinuumProperties:
    photo_ionization_rate_coefficient: pd.DataFrame
    # this is not the rate coefficient but misses Phi I_K
    stimulated_recombination_rate_factor: pd.DataFrame
