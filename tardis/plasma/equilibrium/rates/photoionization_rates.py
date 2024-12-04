import numpy as np
import pandas as pd

from tardis import constants as const
from tardis.transport.montecarlo.estimators.util import (
    integrate_array_by_blocks,
)

C = const.c.cgs.value
H = const.h.cgs.value
K_B = const.k_B.cgs.value


class PhotoionizationRateSolver:
    def __init__(
        self,
        photoionization_cross_sections,
    ):
        self.photoionization_cross_sections = photoionization_cross_sections

        self.photoionization_block_references = np.pad(
            self.photoionization_cross_sections.nu.groupby(level=[0, 1, 2])
            .count()
            .values.cumsum(),
            [1, 0],
        )

        self.photoionization_index = (
            self.photoionization_cross_sections.index.unique()
        )

        self.frequency_i = (
            self.photoionization_cross_sections.groupby(level=[0, 1, 2])
            .first()
            .nu
        )

    def solve(
        self,
        photoionization_rate_estimator,
        stimulated_recombination_rate,
        electron_temperature,
    ):
        """
        Prepares the ionization and recombination coefficients by grouping them for
        ion numbers.

        Parameters
        ----------
        photoionization_rate_estimator : pandas.DataFrame
            The photoionization rate estimator from MCRT.
        stimulated_recombination_rate : pandas.DataFrame
            The stimulated recombination rate from MCRT.
        electron_temperature : u.Quantity
            Electron temperature in each shell.

        Returns
        -------
        photoionization_rate
            Photoionization rate grouped by atomic number and ion number.
        recombination_rate
            Radiative recombination rate grouped by atomic number and ion number.
        """
        nu = self.photoionization_cross_sections["nu"].values
        boltzmann_factor = np.exp(
            -nu[np.newaxis].T / electron_temperature * (H / K_B)
        )

        x_sect = self.photoionization_cross_sections["x_sect"].values
        factor = (
            1 - self.frequency_i / self.photoionization_cross_sections["nu"]
        ).values
        spontaneous_recombination_rate = (
            8 * np.pi * x_sect * factor * nu**3 / C**2
        ) * H
        spontaneous_recombination_rate = (
            spontaneous_recombination_rate[:, np.newaxis] * boltzmann_factor
        )
        spontaneous_recombination_rate = integrate_array_by_blocks(
            spontaneous_recombination_rate,
            nu,
            self.photoionization_block_references,
        )
        spontaneous_recombination_rate = pd.DataFrame(
            spontaneous_recombination_rate, index=self.photoionization_index
        )

        photoionization_rate = photoionization_rate_estimator.groupby(
            level=("atomic_number", "ion_number")
        ).sum()

        recombination_rate = (
            (spontaneous_recombination_rate + stimulated_recombination_rate)
            .groupby(level=["atomic_number", "ion_number"])
            .sum()
        )
        return (
            photoionization_rate,
            recombination_rate,
        )
