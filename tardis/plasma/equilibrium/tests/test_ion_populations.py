from pathlib import Path
from types import SimpleNamespace

import astropy.units as u
import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
from tardisbase.testing.regression_data.regression_data import RegressionData

from tardis.iip_plasma.properties.ion_population import (
    NLTEIonNumberDensity as IIPNLTEIonNumberDensity,
)
from tardis.io.atom_data import AtomData
from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.ion_populations import IonPopulationSolver
from tardis.plasma.equilibrium.rate_matrix import IonRateMatrix
from tardis.plasma.equilibrium.rates import (
    AnalyticPhotoionizationRateSolver,
    CollisionalIonizationRateSolver,
)
from tardis.plasma.radiation_field import (
    DilutePlanckianRadiationField,
)


def hydrogen_population_inputs() -> dict[str, object]:
    """Return small Hydrogen inputs for ion-population solver tests."""
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
    return {
        "radiation_field": radiation_field,
        "thermal_electron_energy_distribution": thermal_electron_energy_distribution,
        "lte_level_population": lte_level_population,
        "estimated_level_population": lte_level_population.copy() * 1.4,
        "lte_ion_population": lte_ion_population,
        "estimated_ion_population": lte_ion_population.copy() * 1.1,
        "partition_function": 1.0,
        "boltzmann_factor": boltzmann_factor,
        "elemental_number_density": elemental_number_density,
    }


def solve_population(
    rate_matrix_solver: IonRateMatrix,
    inputs: dict[str, object],
    charge_conservation: bool,
) -> tuple[pd.DataFrame, pd.Series, IonPopulationSolver]:
    """Solve ion populations and return the solver for matrix inspection."""
    ion_population_solver = IonPopulationSolver(rate_matrix_solver)
    ion_population, electron_density = ion_population_solver.solve(
        inputs["radiation_field"],
        inputs["thermal_electron_energy_distribution"],
        inputs["elemental_number_density"],
        inputs["lte_level_population"],
        inputs["estimated_level_population"],
        inputs["lte_ion_population"],
        inputs["estimated_ion_population"],
        inputs["partition_function"],
        inputs["boltzmann_factor"],
        charge_conservation,
    )
    return ion_population, electron_density, ion_population_solver


def assert_charge_conservation(
    ion_population: pd.DataFrame, electron_density: pd.Series
) -> None:
    """Assert that ion charges reconstruct the electron density."""
    electron_density_from_ions = (
        ion_population
        * ion_population.index.get_level_values("ion_number").to_numpy()[
            :, None
        ]
    ).sum()
    npt.assert_allclose(
        electron_density.to_numpy(),
        electron_density_from_ions.to_numpy(),
        rtol=1e-12,
    )


def assert_hydrogen_matrix_balance(
    ion_population_solver: IonPopulationSolver,
    ion_population: pd.DataFrame,
    elemental_number_density: pd.DataFrame,
) -> None:
    """Assert final Hydrogen matrices satisfy the normalized balance system."""
    for shell in ion_population.columns:
        matrix = ion_population_solver.rates_matrices.loc[1, shell]
        population = ion_population[shell].to_numpy()
        balance = np.array([0.0, elemental_number_density.loc[1, shell]])
        npt.assert_allclose(
            matrix @ population, balance, rtol=1e-12, atol=1e-12
        )
        assert np.isfinite(np.linalg.cond(matrix))


def test_solve(
    rate_matrix_solver: IonRateMatrix, regression_data: RegressionData
) -> None:
    inputs = hydrogen_population_inputs()
    actual_ion_population, actual_electron_density, ion_population_solver = (
        solve_population(rate_matrix_solver, inputs, charge_conservation=False)
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
    npt.assert_allclose(
        actual_ion_population.groupby(level="atomic_number").sum().to_numpy(),
        inputs["elemental_number_density"].to_numpy(),
        rtol=1e-12,
    )
    assert_charge_conservation(actual_ion_population, actual_electron_density)
    assert_hydrogen_matrix_balance(
        ion_population_solver,
        actual_ion_population,
        inputs["elemental_number_density"],
    )


def test_charge_conserving_hydrogen_matches_regression_and_analytic_root(
    rate_matrix_solver: IonRateMatrix, regression_data: RegressionData
) -> None:
    inputs = hydrogen_population_inputs()

    ion_population, electron_density, ion_population_solver = solve_population(
        rate_matrix_solver, inputs, charge_conservation=True
    )

    expected_ion_population = regression_data.sync_dataframe(
        ion_population, key="charge_conserving_ion_population"
    )
    expected_electron_density = regression_data.sync_dataframe(
        electron_density, key="charge_conserving_electron_density"
    )
    pdt.assert_frame_equal(
        ion_population, expected_ion_population, atol=0.0, rtol=1e-12
    )
    pdt.assert_series_equal(
        electron_density, expected_electron_density, atol=0.0, rtol=1e-12
    )

    assert np.all(ion_population.to_numpy() >= 0.0)
    npt.assert_allclose(
        ion_population.groupby(level="atomic_number").sum().to_numpy(),
        inputs["elemental_number_density"].to_numpy(),
        rtol=1e-12,
    )
    assert_charge_conservation(ion_population, electron_density)
    assert_hydrogen_matrix_balance(
        ion_population_solver, ion_population, inputs["elemental_number_density"]
    )

    for shell in ion_population.columns:
        matrix = ion_population_solver.rates_matrices.loc[1, shell]
        photoionization_rate = -matrix[0, 0]
        recombination_coefficient = matrix[0, 1] / electron_density.loc[shell]
        hydrogen_density = inputs["elemental_number_density"].loc[1, shell]
        expected_electron_density = (
            -photoionization_rate
            + np.sqrt(
                photoionization_rate**2
                + 4
                * recombination_coefficient
                * photoionization_rate
                * hydrogen_density
            )
        ) / (2 * recombination_coefficient)
        npt.assert_allclose(
            electron_density.loc[shell],
            expected_electron_density,
            rtol=1e-10,
            atol=0.0,
        )
        npt.assert_allclose(
            ion_population.loc[(1, 1), shell],
            electron_density.loc[shell],
            rtol=1e-12,
            atol=0.0,
        )


def test_charge_conserving_hydrogen_is_seed_independent(
    rate_matrix_solver: IonRateMatrix,
) -> None:
    low_seed_inputs = hydrogen_population_inputs()
    high_seed_inputs = hydrogen_population_inputs()
    low_seed_inputs[
        "thermal_electron_energy_distribution"
    ] = ThermalElectronEnergyDistribution(
        0, np.ones(20) * 10000 * u.K, np.ones(20) * 1.0e2 * u.cm**-3
    )
    high_seed_inputs[
        "thermal_electron_energy_distribution"
    ] = ThermalElectronEnergyDistribution(
        0, np.ones(20) * 10000 * u.K, np.ones(20) * 9.0e4 * u.cm**-3
    )

    low_seed_ions, low_seed_electrons, _ = solve_population(
        rate_matrix_solver, low_seed_inputs, charge_conservation=True
    )
    high_seed_ions, high_seed_electrons, _ = solve_population(
        rate_matrix_solver, high_seed_inputs, charge_conservation=True
    )

    pdt.assert_frame_equal(low_seed_ions, high_seed_ions, rtol=1e-10)
    pdt.assert_series_equal(low_seed_electrons, high_seed_electrons, rtol=1e-10)


def test_charge_conserving_zero_element_density_returns_zero_solution(
    rate_matrix_solver: IonRateMatrix,
) -> None:
    inputs = hydrogen_population_inputs()
    inputs["elemental_number_density"] = inputs[
        "elemental_number_density"
    ].copy()
    inputs["elemental_number_density"].loc[1] = 0.0
    inputs[
        "thermal_electron_energy_distribution"
    ] = ThermalElectronEnergyDistribution(
        0, np.ones(20) * 10000 * u.K, np.ones(20) * 1.0e2 * u.cm**-3
    )

    ion_population, electron_density, _ = solve_population(
        rate_matrix_solver, inputs, charge_conservation=True
    )

    npt.assert_allclose(electron_density.to_numpy(), 0.0)
    npt.assert_allclose(ion_population.to_numpy(), 0.0)


def h_non_h_population_inputs(
    tardis_regression_path: Path,
) -> dict[str, object]:
    """Return H plus one non-H element for real ionization solver tests."""
    columns = pd.Index(["inner", "outer"], name="shell")
    atom_data = AtomData.from_hdf(
        tardis_regression_path
        / "atom_data"
        / "nlte_atom_data"
        / "TestNLTE_He_Ti.h5"
    )
    photoionization_cross_sections = (
        atom_data.photoionization_data.loc[
            [(1, 0, 0), (2, 0, 0), (2, 1, 0)]
        ]
        .sort_values(["atomic_number", "ion_number", "level_number", "nu"])
    )
    level_index = photoionization_cross_sections.index.unique()
    lte_level_population = pd.DataFrame(
        np.ones((len(level_index), len(columns))) * 1.0e5,
        index=level_index,
        columns=columns,
    )
    lte_ion_population = pd.DataFrame(
        np.ones((3, len(columns))) * 1.0e5,
        index=pd.MultiIndex.from_tuples(
            [(1, 1), (2, 1), (2, 2)],
            names=["atomic_number", "ion_number"],
        ),
        columns=columns,
    )
    elemental_number_density = pd.DataFrame(
        [[1.0e5, 2.0e5], [2.0e4, 3.0e4]],
        index=pd.Index([1, 2], name="atomic_number"),
        columns=columns,
    )
    return {
        "radiation_field": DilutePlanckianRadiationField(
            np.ones(len(columns)) * 10000 * u.K,
            dilution_factor=np.ones(len(columns)) * 0.5,
        ),
        "thermal_electron_energy_distribution": ThermalElectronEnergyDistribution(
            0,
            np.ones(len(columns)) * 10000 * u.K,
            np.array([1.0e4, 2.0e4]) * u.cm**-3,
        ),
        "lte_level_population": lte_level_population,
        "estimated_level_population": lte_level_population.copy(),
        "lte_ion_population": lte_ion_population,
        "estimated_ion_population": lte_ion_population.copy(),
        "partition_function": 1.0,
        "boltzmann_factor": pd.DataFrame(
            np.ones_like(lte_level_population),
            index=level_index,
            columns=columns,
        ),
        "elemental_number_density": elemental_number_density,
        "rate_matrix_solver": IonRateMatrix(
            AnalyticPhotoionizationRateSolver(photoionization_cross_sections),
            CollisionalIonizationRateSolver(photoionization_cross_sections),
        ),
    }


def test_charge_conserving_multi_element_solution_matches_regression(
    regression_data: RegressionData, tardis_regression_path: Path
) -> None:
    inputs = h_non_h_population_inputs(tardis_regression_path)

    ion_population, electron_density, _ = solve_population(
        inputs["rate_matrix_solver"], inputs, charge_conservation=True
    )

    expected_ion_population = regression_data.sync_dataframe(
        ion_population, key="h_non_h_charge_conserving_ion_population"
    )
    expected_electron_density = regression_data.sync_dataframe(
        electron_density, key="h_non_h_charge_conserving_electron_density"
    )
    pdt.assert_frame_equal(
        ion_population, expected_ion_population, atol=0.0, rtol=1e-12
    )
    pdt.assert_series_equal(
        electron_density, expected_electron_density, atol=0.0, rtol=1e-12
    )
    pdt.assert_index_equal(
        ion_population.loc[2].index,
        pd.Index(range(3), name="ion_number"),
    )
    pdt.assert_index_equal(
        ion_population.columns,
        inputs["elemental_number_density"].columns,
    )
    npt.assert_allclose(
        ion_population.groupby(level="atomic_number").sum(),
        inputs["elemental_number_density"],
        rtol=1e-12,
    )
    assert_charge_conservation(ion_population, electron_density)


def test_charge_conserving_hydrogen_matches_iip_nlte_solver(
    rate_matrix_solver: IonRateMatrix,
) -> None:
    inputs = hydrogen_population_inputs()
    ion_population, electron_density, ion_population_solver = solve_population(
        rate_matrix_solver, inputs, charge_conservation=True
    )
    shell = ion_population.columns[0]
    matrix = ion_population_solver.rates_matrices.loc[1, shell]
    level_index = pd.MultiIndex.from_tuples(
        [(1, 0, 0)], names=["atomic_number", "ion_number", "level_number"]
    )
    ion_index = pd.MultiIndex.from_tuples(
        [(1, 1)], names=["atomic_number", "ion_number"]
    )
    columns = pd.Index([0], name="shell")
    hydrogen_density = inputs["elemental_number_density"][[shell]].copy()
    hydrogen_density.columns = columns
    gamma = pd.DataFrame([[-matrix[0, 0]]], index=level_index, columns=columns)
    alpha_sp = pd.DataFrame(
        [[matrix[0, 1] / electron_density.loc[shell]]],
        index=level_index,
        columns=columns,
    )
    zero_level_rate = pd.DataFrame([[0.0]], index=level_index, columns=columns)
    iip_ion_population, iip_electron_density = IIPNLTEIonNumberDensity(
        SimpleNamespace(
            previous_ion_number_density=None,
            previous_electron_densities=None,
            nlte_species=[(1, 0)],
        )
    ).calculate(
        pd.DataFrame([[1.0]], index=ion_index, columns=columns),
        zero_level_rate,
        alpha_sp,
        gamma,
        zero_level_rate,
        zero_level_rate,
        hydrogen_density,
        pd.DataFrame([[1.0]], index=level_index, columns=columns),
    )

    actual_ion_population = ion_population[[shell]].copy()
    actual_ion_population.columns = columns
    actual_electron_density = pd.Series(
        [electron_density.loc[shell]], index=columns
    )
    pdt.assert_frame_equal(actual_ion_population, iip_ion_population, rtol=1e-5)
    pdt.assert_series_equal(actual_electron_density, iip_electron_density, rtol=1e-5)
