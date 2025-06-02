import numpy as np
import pandas as pd
import pytest

from tardis.plasma.equilibrium.ion_populations import IonPopulationSolver


def test_solve(rate_matrix_solver, ions_dataframe):
    ion_population_solver = IonPopulationSolver(
        rate_matrix_solver, ions_dataframe
    )

    radiation_field = None
    thermal_electron_energy_distribution = None
    lte_level_population = pd.DataFrame(
        {"population": [1.0, 2.0]},
        index=pd.MultiIndex.from_tuples(
            [(1, 0), (2, 1)], names=["atomic_number", "ion_number"]
        ),
    )
    level_population = lte_level_population.copy()
    lte_ion_population = lte_level_population.copy()
    ion_population = lte_level_population.copy()
    charge_conservation = None

    result = ion_population_solver.solve(
        radiation_field,
        thermal_electron_energy_distribution,
        lte_level_population,
        level_population,
        lte_ion_population,
        ion_population,
        charge_conservation,
    )

    assert isinstance(result, pd.DataFrame)
    assert not result.isnull().values.any()
    assert (result >= 0).all().all()
