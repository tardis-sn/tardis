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
from tardis.plasma.equilibrium.ion_populations import (
    FixedElectronDensityIonPopulationSolver,
    IonPopulationSolver,
)
from tardis.plasma.equilibrium.rate_matrix import AnalyticIonRateMatrix
from tardis.plasma.equilibrium.rates import (
    AnalyticPhotoionizationRateSolver,
    CollisionalIonizationRateSolver,
)
from tardis.plasma.equilibrium.rates.util import (
    align_ion_population_to_level_population,
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
    level_to_continuum_saha_factor = lte_level_population / (
        align_ion_population_to_level_population(
            lte_ion_population, lte_level_population
        ).to_numpy()
        * thermal_electron_energy_distribution.number_density.value
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
        "level_to_continuum_saha_factor": level_to_continuum_saha_factor,
        "elemental_number_density": elemental_number_density,
    }


def solve_population(
    rate_matrix_solver: AnalyticIonRateMatrix,
    inputs: dict[str, object],
    charge_conservation: bool,
) -> tuple[pd.DataFrame, pd.Series, object]:
    """Solve ion populations and return the solver for matrix inspection."""
    ion_population_solver = (
        IonPopulationSolver(rate_matrix_solver)
        if charge_conservation
        else FixedElectronDensityIonPopulationSolver(rate_matrix_solver)
    )
    solve_kwargs = (
        {
            "level_to_continuum_saha_factor": inputs[
                "level_to_continuum_saha_factor"
            ]
        }
        if charge_conservation
        else {}
    )
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
        **solve_kwargs,
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
    ion_population_solver: object,
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
    rate_matrix_solver: AnalyticIonRateMatrix, regression_data: RegressionData
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


def test_charge_conserving_hydrogen_matches_analytic_root(
    rate_matrix_solver: AnalyticIonRateMatrix,
) -> None:
    inputs = hydrogen_population_inputs()

    ion_population, electron_density, ion_population_solver = solve_population(
        rate_matrix_solver, inputs, charge_conservation=True
    )

    for shell in ion_population.columns:
        matrix = ion_population_solver.rates_matrices.loc[1, shell]
        ionization_rate = -matrix[0, 0]
        recombination_rate = matrix[0, 1]
        hydrogen_density = inputs["elemental_number_density"].loc[1, shell]
        expected_electron_density = (
            hydrogen_density
            * ionization_rate
            / (ionization_rate + recombination_rate)
        )
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


def test_charge_conserving_hydrogen_is_seed_independent_from_near_neutral_density(
    rate_matrix_solver: AnalyticIonRateMatrix,
) -> None:
    low_seed_inputs = hydrogen_population_inputs()
    high_seed_inputs = hydrogen_population_inputs()
    low_seed_inputs[
        "thermal_electron_energy_distribution"
    ] = ThermalElectronEnergyDistribution(
        0, np.ones(20) * 10000 * u.K, np.ones(20) * 1.0e-5 * u.cm**-3
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
    thermal_electron_energy_distribution = ThermalElectronEnergyDistribution(
        0,
        np.ones(len(columns)) * 10000 * u.K,
        np.array([1.0e4, 2.0e4]) * u.cm**-3,
    )
    level_to_continuum_saha_factor = lte_level_population / (
        align_ion_population_to_level_population(
            lte_ion_population, lte_level_population
        ).to_numpy()
        * thermal_electron_energy_distribution.number_density.value
    )
    return {
        "radiation_field": DilutePlanckianRadiationField(
            np.ones(len(columns)) * 10000 * u.K,
            dilution_factor=np.ones(len(columns)) * 0.5,
        ),
        "thermal_electron_energy_distribution": thermal_electron_energy_distribution,
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
        "level_to_continuum_saha_factor": level_to_continuum_saha_factor,
        "elemental_number_density": elemental_number_density,
        "rate_matrix_solver": AnalyticIonRateMatrix(
            AnalyticPhotoionizationRateSolver(photoionization_cross_sections),
            CollisionalIonizationRateSolver(photoionization_cross_sections),
        ),
    }


def rate_dataframe_to_level_dataframe(
    rate_dataframe: pd.DataFrame,
) -> pd.DataFrame:
    """Return a level-indexed rate DataFrame from an ion-transition rate."""
    level_rate_dataframe = rate_dataframe.reset_index().set_index(
        ["atomic_number", "ion_number", "level_number_source"]
    )[rate_dataframe.columns]
    level_rate_dataframe.index = level_rate_dataframe.index.set_names(
        ["atomic_number", "ion_number", "level_number"]
    )
    return level_rate_dataframe


def calculate_iip_rate_coefficients(
    rate_matrix_solver: AnalyticIonRateMatrix,
    thermal_electron_energy_distribution: ThermalElectronEnergyDistribution,
    radiation_field: DilutePlanckianRadiationField,
    lte_level_population: pd.DataFrame,
    lte_ion_population: pd.DataFrame,
    estimated_level_population: pd.DataFrame,
    estimated_ion_population: pd.DataFrame,
    partition_function: pd.DataFrame | float,
    boltzmann_factor: pd.DataFrame,
    level_to_continuum_saha_factor: pd.DataFrame,
    electron_density: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Calculate IIP coefficients from the production rate solvers.

    The equilibrium solver returns rates after the electron-density factors in
    Lucy (2003) have been applied. IIP accepts the corresponding coefficients,
    so these rates are converted at the converged electron density.
    """
    photoionization_level_index = (
        rate_matrix_solver.radiative_ionization_rate_solver.photoionization_cross_sections.index.unique()
    )
    lte_level_population = lte_level_population.loc[
        photoionization_level_index
    ]
    estimated_level_population = estimated_level_population.loc[
        photoionization_level_index
    ]
    boltzmann_factor = boltzmann_factor.loc[photoionization_level_index]
    level_to_continuum_saha_factor = level_to_continuum_saha_factor.loc[
        photoionization_level_index
    ]
    final_electron_distribution = ThermalElectronEnergyDistribution(
        thermal_electron_energy_distribution.energy,
        thermal_electron_energy_distribution.temperature,
        electron_density.to_numpy() * u.cm**-3,
    )
    photoionization_rates, spontaneous_recombination_rates = (
        rate_matrix_solver.radiative_ionization_rate_solver.solve(
            radiation_field,
            final_electron_distribution,
            lte_level_population,
            estimated_level_population,
            lte_ion_population,
            estimated_ion_population,
            partition_function,
            boltzmann_factor,
            level_to_continuum_saha_factor,
        )
    )
    collisional_ionization_rates, collisional_recombination_rates = (
        rate_matrix_solver.collisional_ionization_rate_solver.solve(
            final_electron_distribution,
            level_to_continuum_saha_factor,
            partition_function,
            boltzmann_factor,
        )
    )
    level_population_fraction = boltzmann_factor.groupby(
        level=["atomic_number", "ion_number"]
    ).transform(lambda values: values / values.sum())
    gamma = rate_dataframe_to_level_dataframe(photoionization_rates).divide(
        level_population_fraction, axis="index"
    )
    alpha_sp = rate_dataframe_to_level_dataframe(
        spontaneous_recombination_rates
    ).divide(electron_density, axis="columns")
    coll_ion_coeff = rate_dataframe_to_level_dataframe(
        collisional_ionization_rates
    ).divide(electron_density, axis="columns").divide(
        level_population_fraction, axis="index"
    )
    coll_recomb_coeff = rate_dataframe_to_level_dataframe(
        collisional_recombination_rates
    ).divide(electron_density**2, axis="columns")
    return gamma, alpha_sp, coll_ion_coeff, coll_recomb_coeff


def add_missing_iip_ion_levels(
    boltzmann_factor: pd.DataFrame, ion_index: pd.MultiIndex
) -> pd.DataFrame:
    """Add ground levels for ion stages absent from the level data."""
    present_ions = set(
        boltzmann_factor.index.droplevel("level_number").unique().tolist()
    )
    missing_ions = [ion for ion in ion_index if ion not in present_ions]
    if not missing_ions:
        return boltzmann_factor

    missing_level_index = pd.MultiIndex.from_tuples(
        [(atomic_number, ion_number, 0) for atomic_number, ion_number in missing_ions],
        names=boltzmann_factor.index.names,
    )
    return pd.concat(
        [
            boltzmann_factor,
            pd.DataFrame(
                1.0,
                index=missing_level_index,
                columns=boltzmann_factor.columns,
            ),
        ]
    )


def test_charge_conserving_multi_element_solution_uses_real_atomic_data(
    tardis_regression_path: Path,
) -> None:
    inputs = h_non_h_population_inputs(tardis_regression_path)

    ion_population, electron_density, _ = solve_population(
        inputs["rate_matrix_solver"], inputs, charge_conservation=True
    )

    columns = pd.RangeIndex(len(ion_population.columns))
    ion_index = ion_population.index
    gamma, alpha_sp, coll_ion_coeff, coll_recomb_coeff = (
        calculate_iip_rate_coefficients(
            inputs["rate_matrix_solver"],
            inputs["thermal_electron_energy_distribution"],
            inputs["radiation_field"],
            inputs["lte_level_population"],
            inputs["lte_ion_population"],
            inputs["estimated_level_population"],
            inputs["estimated_ion_population"],
            inputs["partition_function"],
            inputs["boltzmann_factor"],
            inputs["level_to_continuum_saha_factor"],
            electron_density,
        )
    )
    iip_level_boltzmann_factor = add_missing_iip_ion_levels(
        inputs["boltzmann_factor"], ion_index
    ).set_axis(columns, axis="columns")
    zero_level_rate = pd.DataFrame(0.0, index=gamma.index, columns=columns)
    iip_ion_population, iip_electron_density = IIPNLTEIonNumberDensity(
        SimpleNamespace(
            previous_ion_number_density=None,
            previous_electron_densities=None,
            nlte_species=[
                (atomic_number, ion_number)
                for atomic_number, ion_number in gamma.index.droplevel(
                    "level_number"
                ).unique()
            ],
        )
    ).calculate(
        pd.DataFrame(
            1.0,
            index=ion_index[ion_index.get_level_values("ion_number") > 0],
            columns=columns,
        ),
        zero_level_rate,
        alpha_sp.set_axis(columns, axis="columns"),
        gamma.set_axis(columns, axis="columns"),
        coll_ion_coeff.set_axis(columns, axis="columns"),
        coll_recomb_coeff.set_axis(columns, axis="columns"),
        inputs["elemental_number_density"].set_axis(columns, axis="columns"),
        iip_level_boltzmann_factor,
    )

    actual_ion_population = ion_population.set_axis(columns, axis="columns")
    actual_electron_density = electron_density.set_axis(columns)
    actual_ion_fraction = actual_ion_population.groupby(
        level="atomic_number"
    ).transform(lambda population: population / population.sum())
    iip_ion_fraction = iip_ion_population.groupby(
        level="atomic_number"
    ).transform(lambda population: population / population.sum())
    pdt.assert_frame_equal(
        actual_ion_fraction, iip_ion_fraction, rtol=1e-5, atol=1e-12
    )
    pdt.assert_series_equal(
        actual_electron_density, iip_electron_density, rtol=1e-5, atol=1e-20
    )


def test_charge_conserving_hydrogen_matches_iip_nlte_solver(
    rate_matrix_solver: AnalyticIonRateMatrix,
) -> None:
    inputs = hydrogen_population_inputs()
    ion_population, electron_density, _ = solve_population(
        rate_matrix_solver, inputs, charge_conservation=True
    )
    shell = ion_population.columns[0]
    phi_index = pd.MultiIndex.from_tuples(
        [(1, 1)], names=["atomic_number", "ion_number"]
    )
    columns = pd.Index([0])
    hydrogen_density = inputs["elemental_number_density"][[shell]].copy()
    hydrogen_density.columns = columns
    gamma, alpha_sp, coll_ion_coeff, coll_recomb_coeff = (
        calculate_iip_rate_coefficients(
            rate_matrix_solver,
            inputs["thermal_electron_energy_distribution"],
            inputs["radiation_field"],
            inputs["lte_level_population"],
            inputs["lte_ion_population"],
            inputs["estimated_level_population"],
            inputs["estimated_ion_population"],
            inputs["partition_function"],
            inputs["boltzmann_factor"],
            inputs["level_to_continuum_saha_factor"],
            electron_density,
        )
    )
    gamma = gamma[[shell]].set_axis(columns, axis="columns")
    alpha_sp = alpha_sp[[shell]].set_axis(columns, axis="columns")
    coll_ion_coeff = coll_ion_coeff[[shell]].set_axis(columns, axis="columns")
    coll_recomb_coeff = coll_recomb_coeff[[shell]].set_axis(
        columns, axis="columns"
    )
    zero_level_rate = pd.DataFrame(0.0, index=gamma.index, columns=columns)
    iip_ion_population, iip_electron_density = IIPNLTEIonNumberDensity(
        SimpleNamespace(
            previous_ion_number_density=None,
            previous_electron_densities=None,
            nlte_species=[(1, 0)],
        )
    ).calculate(
        pd.DataFrame([[1.0]], index=phi_index, columns=columns),
        zero_level_rate,
        alpha_sp,
        gamma,
        coll_ion_coeff,
        coll_recomb_coeff,
        hydrogen_density,
        inputs["boltzmann_factor"][[shell]].set_axis(columns, axis="columns"),
    )

    actual_ion_population = ion_population[[shell]].copy()
    actual_ion_population.columns = columns
    actual_electron_density = pd.Series(
        [electron_density.loc[shell]], index=columns
    )
    actual_ion_fraction = actual_ion_population / actual_ion_population.sum()
    iip_ion_fraction = iip_ion_population / iip_ion_population.sum()
    pdt.assert_frame_equal(
        actual_ion_fraction, iip_ion_fraction, rtol=1e-5, atol=1e-12
    )
    pdt.assert_series_equal(
        actual_electron_density, iip_electron_density, rtol=1e-5, atol=1e-20
    )
