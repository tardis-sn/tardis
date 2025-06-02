import astropy.units as u
import numpy as np
import pandas as pd

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.ion_populations import IonPopulationSolver
from tardis.plasma.radiation_field import (
    DilutePlanckianRadiationField,
)


def test_solve(rate_matrix_solver):
    ion_population_solver = IonPopulationSolver(rate_matrix_solver)

    radiation_field = DilutePlanckianRadiationField(
        np.ones(20) * 10000 * u.K, dilution_factor=np.ones(20) * 0.5
    )
    thermal_electron_energy_distribution = ThermalElectronEnergyDistribution(
        0, np.ones(20) * 10000 * u.K, np.ones(20) * 2e9 * u.cm**-3
    )
    lte_level_population = pd.DataFrame(
        {"population": [1.0, 2.0]},
        index=pd.MultiIndex.from_tuples(
            [(1, 0), (2, 1)], names=["atomic_number", "ion_number"]
        ),
    )
    level_population = lte_level_population.copy()
    lte_ion_population = lte_level_population.copy()
    ion_population = lte_level_population.copy()
    charge_conservation = False

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
