import astropy.units as u
import numpy as np
import pandas as pd
import pandas.testing as pdt

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.ion_populations import IonPopulationSolver
from tardis.plasma.radiation_field import (
    DilutePlanckianRadiationField,
)


def test_solve(rate_matrix_solver, regression_data):
    ion_population_solver = IonPopulationSolver(rate_matrix_solver)

    radiation_field = DilutePlanckianRadiationField(
        np.ones(20) * 10000 * u.K, dilution_factor=np.ones(20) * 0.5
    )
    thermal_electron_energy_distribution = ThermalElectronEnergyDistribution(
        0, np.ones(20) * 10000 * u.K, np.ones(20) * 2e9 * u.cm**-3
    )
    lte_level_population = pd.DataFrame(
        data=np.ones((2, 20)) * 1e5,
        index=pd.MultiIndex.from_tuples(
            [(1, 0, 0), (1, 0, 1)],
            names=["atomic_number", "ion_number", "level_number"],
        ),
    )

    lte_ion_population = pd.DataFrame(
        data=np.vstack([np.ones(20) * 1e5, np.ones(20) * 1e10]),
        index=pd.MultiIndex.from_tuples(
            [(1, 0), (1, 1)],
            names=["atomic_number", "ion_number"],
        ),
    )

    level_population = lte_level_population.copy() * 1.4
    ion_population = lte_ion_population.copy() * 3.0
    # needs to be true for the solver to work
    charge_conservation = True

    actual_ion_population, actual_electron_density = (
        ion_population_solver.solve(
            radiation_field,
            thermal_electron_energy_distribution,
            lte_level_population,
            level_population,
            lte_ion_population,
            ion_population,
            charge_conservation,
        )
    )

    expected_ion_population = regression_data.sync_dataframe(
        actual_ion_population
    )
    expected_electron_density = regression_data.sync_dataframe(
        actual_electron_density
    )

    pdt.assert_frame_equal(actual_ion_population, expected_ion_population)
    pdt.assert_series_equal(actual_electron_density, expected_electron_density)
