from pathlib import Path

import astropy.units as u
import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import pytest
from tardisbase.testing.regression_data.regression_data import RegressionData

from tardis.io.atom_data import AtomData
from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.ion_number_densities import (
    FixedElectronDensityIonNumberDensitySolver,
    IonNumberDensitySolver,
)
from tardis.plasma.equilibrium.rate_matrix import AnalyticIonRateMatrix
from tardis.plasma.equilibrium.rates import (
    AnalyticPhotoionizationRateSolver,
    CollisionalIonizationRateSolver,
)
from tardis.plasma.equilibrium.rates.util import (
    reindex_ion_number_density_to_level_number_density,
)
from tardis.plasma.radiation_field import DilutePlanckianRadiationField


@pytest.fixture
def hydrogen_number_density_inputs() -> dict:
    """Return small Hydrogen inputs for ion-number density solver tests."""
    radiation_field = DilutePlanckianRadiationField(
        np.ones(20) * 10000 * u.K, dilution_factor=np.ones(20) * 0.5
    )
    thermal_electron_energy_distribution = ThermalElectronEnergyDistribution(
        0, np.ones(20) * 10000 * u.K, np.ones(20) * 2e9 * u.cm**-3
    )
    lte_level_number_density = pd.DataFrame(
        data=np.vstack(
            [np.ones(20) * 1e5, np.ones(20) * 1e-1, np.ones(20) * 1e10]
        ),
        index=pd.MultiIndex.from_tuples(
            [(1, 0, 0), (1, 0, 1), (1, 1, 0)],
            names=["atomic_number", "ion_number", "level_number"],
        ),
    )
    lte_ion_number_density = pd.DataFrame(
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
    level_to_continuum_saha_factor = lte_level_number_density / (
        reindex_ion_number_density_to_level_number_density(
            lte_ion_number_density, lte_level_number_density
        ).to_numpy()
        * thermal_electron_energy_distribution.number_density.value
    )
    return {
        "radiation_field": radiation_field,
        "thermal_electron_energy_distribution": thermal_electron_energy_distribution,
        "lte_level_number_density": lte_level_number_density,
        "estimated_level_number_density": lte_level_number_density.copy() * 1.4,
        "lte_ion_number_density": lte_ion_number_density,
        "estimated_ion_number_density": lte_ion_number_density.copy() * 1.1,
        "partition_function": pd.DataFrame(
            1.0,
            index=lte_level_number_density.index.droplevel("level_number").unique(),
            columns=lte_ion_number_density.columns,
        ),
        "boltzmann_factor": boltzmann_factor,
        "level_to_continuum_saha_factor": level_to_continuum_saha_factor,
        "elemental_number_density": elemental_number_density,
    }


def solve_number_density(
    rate_matrix_solver: AnalyticIonRateMatrix,
    inputs: dict,
    charge_conservation: bool,
) -> tuple[pd.DataFrame, pd.Series, object]:
    """Solve ion number densities and return the solver for matrix inspection."""
    ion_number_density_solver = (
        IonNumberDensitySolver(rate_matrix_solver)
        if charge_conservation
        else FixedElectronDensityIonNumberDensitySolver(rate_matrix_solver)
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
    ion_number_density, electron_density = ion_number_density_solver.solve(
        inputs["radiation_field"],
        inputs["thermal_electron_energy_distribution"],
        inputs["elemental_number_density"],
        inputs["lte_level_number_density"],
        inputs["estimated_level_number_density"],
        inputs["lte_ion_number_density"],
        inputs["estimated_ion_number_density"],
        inputs["partition_function"],
        inputs["boltzmann_factor"],
        **solve_kwargs,
    )
    return ion_number_density, electron_density, ion_number_density_solver


def assert_charge_conservation(
    ion_number_density: pd.DataFrame, electron_density: pd.Series
) -> None:
    """Assert that ion charges reconstruct the electron density."""
    electron_density_from_ions = (
        ion_number_density
        * ion_number_density.index.get_level_values("ion_number").to_numpy()[
            :, None
        ]
    ).sum()
    npt.assert_allclose(
        electron_density.to_numpy(),
        electron_density_from_ions.to_numpy(),
        rtol=1e-12,
    )


def assert_hydrogen_matrix_balance(
    ion_number_density_solver: object,
    ion_number_density: pd.DataFrame,
    elemental_number_density: pd.DataFrame,
) -> None:
    """Assert final Hydrogen matrices satisfy the normalized balance system."""
    for shell in ion_number_density.columns:
        matrix = ion_number_density_solver.rates_matrices.loc[1, shell]
        number_density = ion_number_density[shell].to_numpy()
        balance = np.array([0.0, elemental_number_density.loc[1, shell]])
        npt.assert_allclose(
            matrix @ number_density, balance, rtol=1e-12, atol=1e-12
        )
        assert np.isfinite(np.linalg.cond(matrix))


def test_solve(
    rate_matrix_solver: AnalyticIonRateMatrix,
    regression_data: RegressionData,
    hydrogen_number_density_inputs: dict,
) -> None:
    actual_ion_number_density, actual_electron_density, ion_number_density_solver = (
        solve_number_density(
            rate_matrix_solver,
            hydrogen_number_density_inputs,
            charge_conservation=False,
        )
    )

    expected_ion_number_density = regression_data.sync_dataframe(
        actual_ion_number_density, key="ion_number_density"
    )
    expected_electron_density = regression_data.sync_dataframe(
        actual_electron_density, key="electron_density"
    )

    pdt.assert_frame_equal(
        actual_ion_number_density, expected_ion_number_density, atol=0, rtol=1e-15
    )
    pdt.assert_series_equal(
        actual_electron_density, expected_electron_density, atol=0, rtol=1e-15
    )
    assert np.all(actual_ion_number_density.to_numpy() >= 0.0)
    npt.assert_allclose(
        actual_ion_number_density.groupby(level="atomic_number").sum().to_numpy(),
        hydrogen_number_density_inputs["elemental_number_density"].to_numpy(),
        rtol=1e-12,
    )
    assert_charge_conservation(actual_ion_number_density, actual_electron_density)
    assert_hydrogen_matrix_balance(
        ion_number_density_solver,
        actual_ion_number_density,
        hydrogen_number_density_inputs["elemental_number_density"],
    )


def test_charge_conserving_hydrogen_matches_analytic_root(
    rate_matrix_solver: AnalyticIonRateMatrix,
    hydrogen_number_density_inputs: dict,
) -> None:
    ion_number_density, electron_density, ion_number_density_solver = solve_number_density(
        rate_matrix_solver, hydrogen_number_density_inputs, charge_conservation=True
    )

    for shell in ion_number_density.columns:
        matrix = ion_number_density_solver.rates_matrices.loc[1, shell]
        ionization_rate = -matrix[0, 0]
        recombination_rate = matrix[0, 1]
        hydrogen_density = hydrogen_number_density_inputs[
            "elemental_number_density"
        ].loc[1, shell]
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
            ion_number_density.loc[(1, 1), shell],
            electron_density.loc[shell],
            rtol=1e-12,
            atol=0.0,
        )


def test_charge_conserving_hydrogen_is_seed_independent_from_near_neutral_density(
    rate_matrix_solver: AnalyticIonRateMatrix,
    hydrogen_number_density_inputs: dict,
) -> None:
    low_seed_inputs = hydrogen_number_density_inputs.copy()
    high_seed_inputs = hydrogen_number_density_inputs.copy()
    low_seed_inputs["thermal_electron_energy_distribution"] = (
        ThermalElectronEnergyDistribution(
            0, np.ones(20) * 10000 * u.K, np.ones(20) * 1.0e-5 * u.cm**-3
        )
    )
    high_seed_inputs["thermal_electron_energy_distribution"] = (
        ThermalElectronEnergyDistribution(
            0, np.ones(20) * 10000 * u.K, np.ones(20) * 9.0e4 * u.cm**-3
        )
    )

    low_seed_ions, low_seed_electrons, _ = solve_number_density(
        rate_matrix_solver, low_seed_inputs, charge_conservation=True
    )
    high_seed_ions, high_seed_electrons, _ = solve_number_density(
        rate_matrix_solver, high_seed_inputs, charge_conservation=True
    )

    pdt.assert_frame_equal(low_seed_ions, high_seed_ions, rtol=1e-10)
    pdt.assert_series_equal(low_seed_electrons, high_seed_electrons, rtol=1e-10)


@pytest.fixture
def h_non_h_number_density_inputs(tardis_regression_path: Path) -> dict:
    """Return H plus one non-H element for real ionization solver tests."""
    columns = pd.Index(["inner", "outer"], name="shell")
    atom_data = AtomData.from_hdf(
        tardis_regression_path
        / "atom_data"
        / "nlte_atom_data"
        / "TestNLTE_He_Ti.h5"
    )
    photoionization_cross_sections = atom_data.photoionization_data.loc[
        [(1, 0, 0), (2, 0, 0), (2, 1, 0)]
    ].sort_values(["atomic_number", "ion_number", "level_number", "nu"])
    level_index = photoionization_cross_sections.index.unique()
    lte_level_number_density = pd.DataFrame(
        np.ones((len(level_index), len(columns))) * 1.0e5,
        index=level_index,
        columns=columns,
    )
    lte_ion_number_density = pd.DataFrame(
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
    level_to_continuum_saha_factor = lte_level_number_density / (
        reindex_ion_number_density_to_level_number_density(
            lte_ion_number_density, lte_level_number_density
        ).to_numpy()
        * thermal_electron_energy_distribution.number_density.value
    )
    return {
        "radiation_field": DilutePlanckianRadiationField(
            np.ones(len(columns)) * 10000 * u.K,
            dilution_factor=np.ones(len(columns)) * 0.5,
        ),
        "thermal_electron_energy_distribution": thermal_electron_energy_distribution,
        "lte_level_number_density": lte_level_number_density,
        "estimated_level_number_density": lte_level_number_density.copy(),
        "lte_ion_number_density": lte_ion_number_density,
        "estimated_ion_number_density": lte_ion_number_density.copy(),
        "partition_function": pd.DataFrame(
            1.0,
            index=level_index.droplevel("level_number").unique(),
            columns=lte_ion_number_density.columns,
        ),
        "boltzmann_factor": pd.DataFrame(
            np.ones_like(lte_level_number_density),
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


def test_charge_conserving_multi_element_solution_uses_real_atomic_data(
    h_non_h_number_density_inputs: dict,
    regression_data: RegressionData,
) -> None:
    inputs = h_non_h_number_density_inputs
    ion_number_density, electron_density, _ = solve_number_density(
        inputs["rate_matrix_solver"], inputs, charge_conservation=True
    )

    columns = pd.RangeIndex(len(ion_number_density.columns))
    actual_ion_number_density = ion_number_density.set_axis(columns, axis="columns")
    actual_ion_fraction = actual_ion_number_density.groupby(
        level="atomic_number"
    ).transform(lambda number_density: number_density / number_density.sum())
    expected_ion_fraction = regression_data.sync_dataframe(
        actual_ion_fraction, key="frame_0"
    )
    pdt.assert_frame_equal(
        actual_ion_fraction,
        expected_ion_fraction,
        rtol=1e-5,
        atol=1e-12,
    )

    actual_electron_density = electron_density.set_axis(columns).to_frame("value")
    expected_electron_density = regression_data.sync_dataframe(
        actual_electron_density, key="series_0"
    )
    pdt.assert_frame_equal(
        actual_electron_density,
        expected_electron_density,
        rtol=1e-5,
        atol=1e-20,
    )


def test_charge_conserving_hydrogen_matches_iip_nlte_solver(
    rate_matrix_solver: AnalyticIonRateMatrix,
    hydrogen_number_density_inputs: dict,
    regression_data: RegressionData,
) -> None:
    ion_number_density, electron_density, _ = solve_number_density(
        rate_matrix_solver, hydrogen_number_density_inputs, charge_conservation=True
    )
    shell = ion_number_density.columns[0]
    columns = pd.Index([0])
    actual_ion_number_density = ion_number_density[[shell]].set_axis(columns, axis="columns")
    actual_ion_fraction = actual_ion_number_density / actual_ion_number_density.sum()
    expected_ion_fraction = regression_data.sync_dataframe(
        actual_ion_fraction, key="frame_0"
    )
    pdt.assert_frame_equal(
        actual_ion_fraction,
        expected_ion_fraction,
        rtol=1e-5,
        atol=1e-12,
    )

    actual_electron_density = pd.Series(
        [electron_density.loc[shell]], index=columns
    ).to_frame("value")
    expected_electron_density = regression_data.sync_dataframe(
        actual_electron_density, key="series_0"
    )
    pdt.assert_frame_equal(
        actual_electron_density,
        expected_electron_density,
        rtol=1e-5,
        atol=1e-20,
    )
