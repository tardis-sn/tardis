from dataclasses import dataclass

import numpy as np
import pandas as pd
from astropy import units as u

import tardis.constants as const
from tardis.io.atom_data import AtomData
from tardis.transport.montecarlo.estimators.radfield_mc_estimators import (
    RadiationFieldMCEstimators,
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
        radfield_mc_estimators: RadiationFieldMCEstimators,
        time_simulation,
        volume,
    ):
        """
        Solve for the continuum properties.

        Parameters
        ----------
        radfield_mc_estimators : RadiationFieldMCEstimators
            The Monte Carlo estimators for the radiation field.
        time_simulation : float
            The simulation time.
        volume : float
            The volume of the cells.

        Returns
        -------
        ContinuumProperties
            The calculated continuum properties.
        """
        photo_ion_norm_factor = (time_simulation * volume * H) ** -1

        photo_ionization_rate_coefficient = bound_free_estimator_array2frame(
            radfield_mc_estimators.photo_ion_estimator,
            self.atom_data.level2continuum_edge_idx,
        )
        photo_ionization_rate_coefficient *= photo_ion_norm_factor

        stimulated_recomb_rate_factor = bound_free_estimator_array2frame(
            radfield_mc_estimators.stim_recomb_estimator,
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
