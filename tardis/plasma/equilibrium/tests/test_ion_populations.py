import astropy.units as u
import numpy as np
import numpy.testing as npt
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
        data=np.vstack(
            [np.ones(20) * 1e5, np.ones(20) * 1e-1, np.ones(20) * 1e10]
        ),
        index=pd.MultiIndex.from_tuples(
            [(1, 0, 0), (1, 0, 1), (1, 1, 0)],
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

    boltzmann_factor = pd.DataFrame(
        data=np.vstack(
            [np.ones(20) * 2.0, np.ones(20) * 0.000011, np.ones(20)]
        ),
        index=pd.MultiIndex.from_tuples(
            [(1, 0, 0), (1, 0, 1), (1, 1, 0)],
            names=["atomic_number", "ion_number", "level_number"],
        ),
    )

    elemental_number_density = pd.DataFrame(
        data=np.vstack([np.ones(20) * 1e5]),
        index=pd.Index([1], name="atomic_number"),
    )

    level_population = lte_level_population.copy() * 1.4
    ion_population = lte_ion_population.copy() * 1.1
    charge_conservation = False

    actual_ion_population, actual_electron_density = (
        ion_population_solver.solve(
            radiation_field,
            thermal_electron_energy_distribution,
            elemental_number_density,
            lte_level_population,
            level_population,
            lte_ion_population,
            ion_population,
            1.0,
            boltzmann_factor,
            charge_conservation,
        )
    )

    expected_ion_population = regression_data.sync_dataframe(
        actual_ion_population, key="ion_population"
    )
    expected_electron_density = regression_data.sync_dataframe(
        actual_electron_density, key="electron_density"
    )

    pdt.assert_frame_equal(
        actual_ion_population, expected_ion_population, atol=0, rtol=1e-15
    )
    pdt.assert_series_equal(
        actual_electron_density, expected_electron_density, atol=0, rtol=1e-15
    )

    assert np.all(actual_ion_population.to_numpy() >= 0.0)
    number_density_from_ions = actual_ion_population.groupby(
        level="atomic_number"
    ).sum()
    pdt.assert_index_equal(
        number_density_from_ions.index,
        elemental_number_density.index,
        check_names=False,
    )
    npt.assert_allclose(
        number_density_from_ions.to_numpy(),
        elemental_number_density.to_numpy(),
        rtol=1e-12,
    )

    electron_density_from_ions = (
        actual_ion_population
        * actual_ion_population.index.get_level_values("ion_number").to_numpy()[
            :, None
        ]
    ).sum()
    npt.assert_allclose(
        actual_electron_density.to_numpy(),
        electron_density_from_ions.to_numpy(),
        rtol=1e-12,
    )

    for shell in actual_ion_population.columns:
        matrix = ion_population_solver.rates_matrices.loc[1, shell]
        population = actual_ion_population[shell].to_numpy()
        balance = np.array([0.0, elemental_number_density.loc[1, shell]])
        npt.assert_allclose(
            matrix @ population, balance, rtol=1e-12, atol=1e-12
        )
        assert np.isfinite(np.linalg.cond(matrix))
