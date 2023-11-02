import pandas as pd
import numpy as np

from tardis.plasma.properties.base import ProcessingPlasmaProperty

__all__ = [
    "NLTEPopulationSolver",
]


class NLTEPopulationSolver(ProcessingPlasmaProperty):
    outputs = ("ion_number_density", "electron_densities")

    def calculate(
        self,
        gamma,
        alpha_sp,
        alpha_stim,
        coll_ion_coeff,
        coll_recomb_coeff,
        partition_function,
        levels,
        level_boltzmann_factor,
        phi,
        rate_matrix_index,
        number_density,
        nlte_excitation_species,
    ):
        initial_electron_densities = number_density.sum(axis=0)

        electron_densities = pd.Series(0.0, index=phi.columns)
        ion_number_density = pd.DataFrame(
            0.0, index=phi.index, columns=phi.columns
        )

        raise NotImplementedError("This is not implemented yet")

        return ion_number_density, electron_densities
