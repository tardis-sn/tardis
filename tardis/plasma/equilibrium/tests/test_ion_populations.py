from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

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
from tardis.plasma.equilibrium.ion_populations import (
    CHARGE_TOLERANCE,
    FixedElectronDensityIonPopulationSolver,
    IonPopulationSolver,
)
from tardis.plasma.equilibrium.rate_matrix import AnalyticIonRateMatrix
from tardis.plasma.equilibrium.rates import (
    AnalyticPhotoionizationRateSolver,
    CollisionalIonizationRateSolver,
)
from tardis.plasma.equilibrium.rates.util import (
    reindex_ion_population_to_level_population,
)
from tardis.plasma.exceptions import PlasmaIonizationError
from tardis.plasma.radiation_field import (
    DilutePlanckianRadiationField,
)


@pytest.fixture
def hydrogen_population_inputs() -> dict:
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
        reindex_ion_population_to_level_population(
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
        "partition_function": pd.DataFrame(
            1.0,
            index=lte_level_population.index.droplevel("level_number").unique(),
            columns=lte_ion_population.columns,
        ),
        "boltzmann_factor": boltzmann_factor,
        "level_to_continuum_saha_factor": level_to_continuum_saha_factor,
        "elemental_number_density": elemental_number_density,
    }


def solve_population(
    rate_matrix_solver: AnalyticIonRateMatrix,
    inputs: dict,
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
    rate_matrix_solver: AnalyticIonRateMatrix,
    regression_data: RegressionData,
    hydrogen_population_inputs: dict,
) -> None:
    actual_ion_population, actual_electron_density, ion_population_solver = (
        solve_population(
            rate_matrix_solver,
            hydrogen_population_inputs,
            charge_conservation=False,
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
    npt.assert_allclose(
        actual_ion_population.groupby(level="atomic_number").sum().to_numpy(),
        hydrogen_population_inputs["elemental_number_density"].to_numpy(),
        rtol=1e-12,
    )
    assert_charge_conservation(actual_ion_population, actual_electron_density)
    assert_hydrogen_matrix_balance(
        ion_population_solver,
        actual_ion_population,
        hydrogen_population_inputs["elemental_number_density"],
    )


def test_charge_conserving_hydrogen_matches_analytic_root(
    rate_matrix_solver: AnalyticIonRateMatrix,
    hydrogen_population_inputs: dict,
) -> None:
    ion_population, electron_density, ion_population_solver = solve_population(
        rate_matrix_solver, hydrogen_population_inputs, charge_conservation=True
    )

    for shell in ion_population.columns:
        matrix = ion_population_solver.rates_matrices.loc[1, shell]
        ionization_rate = -matrix[0, 0]
        recombination_rate = matrix[0, 1]
        hydrogen_density = hydrogen_population_inputs[
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
            ion_population.loc[(1, 1), shell],
            electron_density.loc[shell],
            rtol=1e-12,
            atol=0.0,
        )


def test_charge_conserving_hydrogen_is_seed_independent_from_near_neutral_density(
    rate_matrix_solver: AnalyticIonRateMatrix,
    hydrogen_population_inputs: dict,
) -> None:
    low_seed_inputs = hydrogen_population_inputs.copy()
    high_seed_inputs = hydrogen_population_inputs.copy()
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

    low_seed_ions, low_seed_electrons, _ = solve_population(
        rate_matrix_solver, low_seed_inputs, charge_conservation=True
    )
    high_seed_ions, high_seed_electrons, _ = solve_population(
        rate_matrix_solver, high_seed_inputs, charge_conservation=True
    )

    pdt.assert_frame_equal(low_seed_ions, high_seed_ions, rtol=1e-10)
    pdt.assert_series_equal(low_seed_electrons, high_seed_electrons, rtol=1e-10)


@pytest.mark.parametrize(
    ("prepared_rate_matrix", "rate_matrix_mode"),
    [(False, "direct"), (True, "prepared")],
)
def test_ion_solver_reports_convergence_blocker(
    caplog: pytest.LogCaptureFixture,
    prepared_rate_matrix: bool,
    rate_matrix_mode: str,
) -> None:
    """Report which residual prevents ion-solver convergence.

    Claim
    -----
    Each ion fixed-point iteration identifies whether population or charge
    residuals prevent its active shells from converging.

    Regime
    ------
    Two hydrogen stages in two shells have converged populations. The outer
    shell has a normalized charge residual of ``-0.1``.

    Verification
    ------------
    All population residuals are exactly zero. The local and full outer-shell
    matrices have exact solutions ``[5, 5]`` and ``[6, 4] cm^-3``, respectively,
    producing a normalized final charge residual of ``-0.1``.
    """
    shell_columns = pd.Index(["inner", "outer"], name="shell")
    ion_index = pd.MultiIndex.from_tuples(
        [(1, 0), (1, 1)], names=["atomic_number", "ion_number"]
    )
    rate_matrix_solver = SimpleNamespace(ion_population_index=ion_index)
    if prepared_rate_matrix:
        rate_matrix_solver.prepare = Mock()
    solver = IonPopulationSolver(rate_matrix_solver, max_solver_iterations=1)
    solver._charge_root_diagnostics = {}
    local_matrix = np.array([[-1.0, 1.0], [1.0, 1.0]])
    full_matrix = np.array([[-1.0, 1.5], [1.0, 1.0]])

    def solve_shell_charge(shell_idx: int, *_: object, **__: object) -> float:
        shell = shell_columns[shell_idx]
        solver._charge_root_diagnostics[shell] = (
            0.0,
            7,
            8,
            True,
            "converged",
        )
        return 5.0

    solver.solve_shell_charge = Mock(side_effect=solve_shell_charge)
    solution_population = np.array([[5.0, 6.0], [5.0, 4.0]])
    full_rate_matrices = pd.DataFrame(
        [[local_matrix, full_matrix]],
        index=pd.Index([1], name="atomic_number"),
        columns=shell_columns,
    )

    def solve_charge_balance(
        *_: object, **__: object
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        solver.rates_matrices = full_rate_matrices
        return solution_population, np.array([0.0, -0.1])

    solver.solve_charge_balance = Mock(side_effect=solve_charge_balance)
    estimated_ion_population = pd.DataFrame(
        solution_population,
        index=ion_index,
        columns=shell_columns,
    )
    empty_frame = pd.DataFrame(columns=shell_columns)
    empty_level_frame = pd.DataFrame(
        columns=shell_columns,
        index=pd.MultiIndex.from_tuples(
            [],
            names=["atomic_number", "ion_number", "level_number"],
        ),
    )
    caplog.set_level(
        "DEBUG", logger="tardis.plasma.equilibrium.ion_populations"
    )

    with pytest.raises(
        PlasmaIonizationError,
        match="Ion population solver did not converge",
    ):
        solver._solve_charge_conserving(
            None,
            ThermalElectronEnergyDistribution(
                0.0 * u.erg,
                np.array([6000.0, 6000.0]) * u.K,
                np.array([5.0, 5.0]) / u.cm**3,
            ),
            pd.DataFrame(
                [[10.0, 10.0]],
                index=pd.Index([1], name="atomic_number"),
                columns=shell_columns,
            ),
            empty_frame,
            empty_level_frame,
            empty_frame,
            estimated_ion_population,
            empty_frame,
            empty_frame,
            empty_frame,
            1.0e-10,
        )

    assert (
        "Ion solver iteration 1: 1/2 shells unconverged; population-failing "
        "shells=0, charge-failing shells=1, failures of both=0; maximum "
        "ion-stage relative residual 0.000e+00 at ion stage (1, 0) in shell "
        "inner (tolerance=1.000e-10); maximum normalized charge residual "
        "1.000e-01 in shell outer (tolerance=1.000e-10); blocker=charge"
        in caplog.text
    )
    assert (
        "Ion solver iteration 1 charge-root consistency in shell outer: "
        "local terminal residual=0.000e+00, final residual=-1.000e-01, "
        "absolute difference=1.000e-01; Brent converged=True, iterations=7, "
        "function_calls=8, flag=converged" in caplog.text
    )
    if prepared_rate_matrix:
        assert (
            "Ion solver iteration 1 starting ion rate-matrix preparation"
            in caplog.text
        )
        assert (
            "Ion solver iteration 1 finished ion rate-matrix preparation in"
            in caplog.text
        )
        rate_matrix_solver.prepare.assert_called_once()
    else:
        assert (
            "Ion solver iteration 1 uses direct ion rate-matrix construction"
            in caplog.text
        )
    assert (
        "Ion solver iteration 1 starting 2 shell charge solves "
        f"(rate-matrix mode={rate_matrix_mode})" in caplog.text
    )
    assert (
        "Ion solver iteration 1 finished 2 shell charge solves in"
        in caplog.text
    )
    assert (
        "Ion solver iteration 1 starting full charge-balance reconstruction"
        in caplog.text
    )
    assert (
        "Ion solver iteration 1 finished full charge-balance reconstruction in"
        in caplog.text
    )
    assert "ion-matrix direct-solve consistency" not in caplog.text
    assert not any(
        record.levelname == "INFO" and "Ion solver iteration" in record.message
        for record in caplog.records
    )


def test_shell_charge_records_terminal_root_diagnostic(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Record the residual and metadata at the returned charge root.

    Claim
    -----
    The shell charge solve evaluates its normalized charge residual at the
    electron-density fraction returned by Brent and retains Brent's status.

    Regime
    ------
    A one-shell analytic charge residual has its exact root at an electron-
    density fraction of ``0.25`` and a maximum density of ``40 cm^-3``.

    Verification
    ------------
    The independent linear residual is exactly zero at ``0.25`` and the
    corresponding electron density is ``10 cm^-3``.
    """
    shell_columns = pd.Index(["shell"], name="shell")
    ion_index = pd.MultiIndex.from_tuples(
        [(1, 0), (1, 1)], names=["atomic_number", "ion_number"]
    )
    solver = IonPopulationSolver(
        SimpleNamespace(ion_population_index=ion_index)
    )
    maximum_electron_density = 40.0

    def solve_charge_balance(
        electron_density: npt.NDArray[np.float64],
        *_: object,
        **__: object,
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        residual = electron_density / maximum_electron_density - 0.25
        solver.rates_matrices = pd.DataFrame(
            [[np.array([[-1.0, 1.0], [1.0, 1.0]])]],
            index=pd.Index([1], name="atomic_number"),
            columns=shell_columns,
        )
        return np.zeros((2, 1)), residual

    def solve_root(
        residual: Callable[[float], float],
        lower_bound: float,
        upper_bound: float,
        **kwargs: object,
    ) -> tuple[float, SimpleNamespace]:
        assert lower_bound < 0.25 < upper_bound
        assert kwargs["full_output"] is True
        assert residual(0.25) == 0.0
        return 0.25, SimpleNamespace(
            iterations=4,
            function_calls=6,
            converged=True,
            flag="converged",
        )

    solver.solve_charge_balance = Mock(side_effect=solve_charge_balance)
    monkeypatch.setattr(
        "tardis.plasma.equilibrium.ion_populations.brentq", solve_root
    )
    empty_frame = pd.DataFrame(columns=shell_columns)
    elemental_number_density = pd.DataFrame(
        [[10.0]],
        index=pd.Index([1], name="atomic_number"),
        columns=shell_columns,
    )

    electron_density = solver.solve_shell_charge(
        0,
        maximum_electron_density,
        None,
        ThermalElectronEnergyDistribution(
            0.0 * u.erg,
            np.array([6000.0]) * u.K,
            np.array([1.0]) / u.cm**3,
        ),
        empty_frame,
        empty_frame,
        empty_frame,
        empty_frame,
        empty_frame,
        empty_frame,
        empty_frame,
        elemental_number_density,
    )

    assert electron_density == 10.0
    assert solver._charge_root_diagnostics["shell"] == (
        0.0,
        4,
        6,
        True,
        "converged",
    )


def test_singular_rate_matrix_reports_element_shell_and_density() -> None:
    """A disconnected shell reports its physical and numerical location.

    Claim
    -----
    A singular ion-stage balance identifies the element, shell, trial electron
    density, and rank deficiency that caused the solve to fail.

    Regime
    ------
    One hydrogen shell has a connected normalized two-stage matrix and the
    other contains only its normalization row.

    Verification
    ------------
    The disconnected matrix has exact rank one, independently establishing
    the expected diagnostic values without relying on the production solver.
    """
    shells = pd.Index(["connected", "disconnected"], name="shell")
    connected_matrix = np.array([[-1.0, 1.0], [1.0, 1.0]])
    disconnected_matrix = np.array([[0.0, 0.0], [1.0, 1.0]])
    assert np.linalg.matrix_rank(disconnected_matrix) == 1
    rate_matrix_solver = Mock()
    rate_matrix_solver.prepared_matrix_coefficients = None
    rate_matrix_solver.ion_population_index = pd.MultiIndex.from_tuples(
        [(1, 0), (1, 1)], names=["atomic_number", "ion_number"]
    )
    rate_matrix_solver.solve.return_value = pd.DataFrame(
        [[connected_matrix, disconnected_matrix]],
        index=pd.Index([1], name="atomic_number"),
        columns=shells,
    )
    solver = IonPopulationSolver(rate_matrix_solver)
    electron_distribution = ThermalElectronEnergyDistribution(
        0.0 * u.erg,
        np.array([5000.0, 6000.0]) * u.K,
        np.array([2.0e7, 3.0e7]) / u.cm**3,
    )
    empty_frame = pd.DataFrame(columns=shells)
    elemental_number_density = pd.DataFrame(
        [[1.0e8, 1.0e8]],
        index=pd.Index([1], name="atomic_number"),
        columns=shells,
    )

    with pytest.raises(
        PlasmaIonizationError,
        match=(
            r"atomic number 1, shell disconnected .*"
            r"electron_density=3\.000e\+07 cm\^-3, rank=1/2"
        ),
    ):
        solver.solve_element_populations_at_electron_density(
            np.array([2.0e7, 3.0e7]),
            None,
            electron_distribution,
            empty_frame,
            empty_frame,
            empty_frame,
            empty_frame,
            empty_frame,
            empty_frame,
            empty_frame,
            elemental_number_density,
            maximum_electron_density=np.array([1.0e8, 1.0e8]),
        )


def test_zero_density_element_skips_singular_shell_matrix() -> None:
    """Return zero ions for an absent element without hiding populated shells.

    Claim
    -----
    An element with exactly zero number density has zero population in every
    ion stage and does not require an invertible rate matrix.

    Regime
    ------
    Hydrogen is absent in one shell with a disconnected matrix and has number
    density ``8e7 cm^-3`` in a second shell with equal forward and reverse
    rates.

    Verification
    ------------
    Conservation requires exact zeros in the empty shell. Direct solution of
    the populated 2-by-2 matrix gives equal fractions and therefore ``4e7
    cm^-3`` in each ion stage.
    """
    shells = pd.Index(["empty", "populated"], name="shell")
    disconnected_matrix = np.array([[0.0, 0.0], [1.0, 1.0]])
    connected_matrix = np.array([[-1.0, 1.0], [1.0, 1.0]])
    rate_matrix_solver = Mock()
    rate_matrix_solver.prepared_matrix_coefficients = None
    rate_matrix_solver.ion_population_index = pd.MultiIndex.from_tuples(
        [(1, 0), (1, 1)], names=["atomic_number", "ion_number"]
    )
    rate_matrix_solver.solve.return_value = pd.DataFrame(
        [[disconnected_matrix, connected_matrix]],
        index=pd.Index([1], name="atomic_number"),
        columns=shells,
    )
    solver = IonPopulationSolver(rate_matrix_solver)
    electron_distribution = ThermalElectronEnergyDistribution(
        0.0 * u.erg,
        np.array([5000.0, 6000.0]) * u.K,
        np.array([1.0e-10, 3.0e7]) / u.cm**3,
    )
    empty_frame = pd.DataFrame(columns=shells)
    elemental_number_density = pd.DataFrame(
        [[0.0, 8.0e7]],
        index=pd.Index([1], name="atomic_number"),
        columns=shells,
    )

    ion_population = solver.solve_element_populations_at_electron_density(
        np.array([1.0e-10, 3.0e7]),
        None,
        electron_distribution,
        empty_frame,
        empty_frame,
        empty_frame,
        empty_frame,
        empty_frame,
        empty_frame,
        empty_frame,
        elemental_number_density,
    )

    npt.assert_array_equal(
        ion_population,
        np.array([[0.0, 4.0e7], [0.0, 4.0e7]]),
    )


def test_missing_element_skips_singular_ion_matrix() -> None:
    """Return zero ion populations when an element has no density row.

    Claim
    -----
    An element omitted from the elemental-density index is absent and has zero
    population in every ion stage.

    Regime
    ------
    The ion-rate data contain a disconnected hydrogen matrix, while the model
    density contains only helium.

    Verification
    ------------
    Element conservation requires exact zero hydrogen populations; the helium
    density must not be mistaken for hydrogen through negative indexing.
    """
    shells = pd.Index([0], name="shell")
    disconnected_matrix = np.array([[0.0, 0.0], [1.0, 1.0]])
    rate_matrix_solver = Mock()
    rate_matrix_solver.prepared_matrix_coefficients = None
    rate_matrix_solver.ion_population_index = pd.MultiIndex.from_tuples(
        [(1, 0), (1, 1)], names=["atomic_number", "ion_number"]
    )
    rate_matrix_solver.solve.return_value = pd.DataFrame(
        [[disconnected_matrix]],
        index=pd.Index([1], name="atomic_number"),
        columns=shells,
    )
    solver = IonPopulationSolver(rate_matrix_solver)
    electron_distribution = ThermalElectronEnergyDistribution(
        0.0 * u.erg,
        np.array([6000.0]) * u.K,
        np.array([2.353e-6]) / u.cm**3,
    )
    empty_frame = pd.DataFrame(columns=shells)
    helium_number_density = pd.DataFrame(
        [[1.0e8]],
        index=pd.Index([2], name="atomic_number"),
        columns=shells,
    )

    ion_population = solver.solve_element_populations_at_electron_density(
        np.array([2.353e-6]),
        None,
        electron_distribution,
        empty_frame,
        empty_frame,
        empty_frame,
        empty_frame,
        empty_frame,
        empty_frame,
        empty_frame,
        helium_number_density,
    )

    npt.assert_array_equal(ion_population, np.zeros((2, 1)))


def test_singular_populated_element_freezes_estimated_ion_fractions(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Freeze a populated element when its ion stages are disconnected.

    Claim
    -----
    A singular rate matrix with a valid pre-evaluator population retains those
    ion fractions, conserves the current elemental density, and reports the
    fallback.

    Regime
    ------
    One populated hydrogen shell has a rank-one matrix and estimated neutral
    and ionized populations in a 3:1 ratio.

    Verification
    ------------
    Normalizing the independent 3:1 seed to ``8e7 cm^-3`` gives exact stage
    populations of ``6e7`` and ``2e7 cm^-3``.
    """
    shells = pd.Index(["disconnected"], name="shell")
    disconnected_matrix = np.array([[0.0, 0.0], [1.0, 1.0]])
    ion_index = pd.MultiIndex.from_tuples(
        [(1, 0), (1, 1)], names=["atomic_number", "ion_number"]
    )
    rate_matrix_solver = Mock()
    rate_matrix_solver.prepared_matrix_coefficients = None
    rate_matrix_solver.ion_population_index = ion_index
    rate_matrix_solver.solve.return_value = pd.DataFrame(
        [[disconnected_matrix]],
        index=pd.Index([1], name="atomic_number"),
        columns=shells,
    )
    solver = IonPopulationSolver(rate_matrix_solver)
    electron_distribution = ThermalElectronEnergyDistribution(
        0.0 * u.erg,
        np.array([6000.0]) * u.K,
        np.array([3.0e7]) / u.cm**3,
    )
    empty_frame = pd.DataFrame(columns=shells)
    estimated_ion_population = pd.DataFrame(
        [[6.0e7], [2.0e7]], index=ion_index, columns=shells
    )
    elemental_number_density = pd.DataFrame(
        [[8.0e7]],
        index=pd.Index([1], name="atomic_number"),
        columns=shells,
    )

    ion_population = solver.solve_element_populations_at_electron_density(
        np.array([3.0e7]),
        None,
        electron_distribution,
        empty_frame,
        empty_frame,
        empty_frame,
        estimated_ion_population,
        empty_frame,
        empty_frame,
        empty_frame,
        elemental_number_density,
    )

    npt.assert_array_equal(ion_population[:, 0], np.array([6.0e7, 2.0e7]))
    assert (
        "Freezing ion populations for atomic number 1, shell disconnected"
        in (caplog.text)
    )


def test_negligible_singular_element_without_seed_is_skipped(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Skip an unresolved element whose maximum charge is negligible.

    Claim
    -----
    During charge balance, a singular element with no usable prior population
    may remain at zero only when its maximum possible charge contribution is
    below the charge-convergence tolerance.

    Regime
    ------
    Copper has ``6.44e-21 cm^-3`` total density, 30 disconnected ion stages,
    and a maximum shell electron density of ``6.80e6 cm^-3``.

    Verification
    ------------
    The independently calculated maximum normalized charge contribution is
    ``29 n_Cu / n_e,max = 2.75e-26``, far below ``CHARGE_TOLERANCE``. All Cu
    stages must remain exactly zero and the approximation must be reported.
    """
    atomic_number = 29
    number_of_ion_stages = atomic_number + 1
    shells = pd.Index([39], name="shell")
    disconnected_matrix = np.zeros((number_of_ion_stages, number_of_ion_stages))
    disconnected_matrix[1] = 1.0
    ion_index = pd.MultiIndex.from_product(
        [[atomic_number], range(number_of_ion_stages)],
        names=["atomic_number", "ion_number"],
    )
    rate_matrix_solver = Mock()
    rate_matrix_solver.prepared_matrix_coefficients = None
    rate_matrix_solver.ion_population_index = ion_index
    rate_matrix_solver.solve.return_value = pd.DataFrame(
        [[disconnected_matrix]],
        index=pd.Index([atomic_number], name="atomic_number"),
        columns=shells,
    )
    solver = IonPopulationSolver(rate_matrix_solver)
    electron_distribution = ThermalElectronEnergyDistribution(
        0.0 * u.erg,
        np.array([7937.46]) * u.K,
        np.array([1.510e-9]) / u.cm**3,
    )
    empty_frame = pd.DataFrame(columns=shells)
    elemental_density = 6.4400828601381825e-21
    maximum_electron_density = np.array([6800639.422190048])
    maximum_charge_fraction = (
        atomic_number * elemental_density / maximum_electron_density[0]
    )
    assert maximum_charge_fraction < CHARGE_TOLERANCE
    elemental_number_density = pd.DataFrame(
        [[elemental_density]],
        index=pd.Index([atomic_number], name="atomic_number"),
        columns=shells,
    )

    ion_population, charge_residual = solver.solve_charge_balance(
        np.array([1.510e-9]),
        None,
        electron_distribution,
        empty_frame,
        empty_frame,
        empty_frame,
        empty_frame,
        empty_frame,
        empty_frame,
        empty_frame,
        elemental_number_density,
        maximum_electron_density,
    )

    npt.assert_array_equal(ion_population, np.zeros((number_of_ion_stages, 1)))
    npt.assert_allclose(
        charge_residual, -np.array([1.510e-9]) / maximum_electron_density
    )
    assert "Skipping ion populations for atomic number 29, shell 39" in (
        caplog.text
    )


@pytest.fixture
def h_non_h_population_inputs(
    tardis_regression_path: Path,
) -> dict:
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
        reindex_ion_population_to_level_population(
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
        "partition_function": pd.DataFrame(
            1.0,
            index=level_index.droplevel("level_number").unique(),
            columns=lte_ion_population.columns,
        ),
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


def test_charge_conserving_multi_element_solution_uses_real_atomic_data(
    h_non_h_population_inputs: dict,
    regression_data: RegressionData,
) -> None:
    inputs = h_non_h_population_inputs
    ion_population, electron_density, _ = solve_population(
        inputs["rate_matrix_solver"], inputs, charge_conservation=True
    )

    columns = pd.RangeIndex(len(ion_population.columns))
    actual_ion_population = ion_population.set_axis(columns, axis="columns")
    actual_ion_fraction = actual_ion_population.groupby(
        level="atomic_number"
    ).transform(lambda population: population / population.sum())
    expected_ion_fraction = regression_data.sync_dataframe(
        actual_ion_fraction, key="frame_0"
    )
    pdt.assert_frame_equal(
        actual_ion_fraction,
        expected_ion_fraction,
        rtol=1e-5,
        atol=1e-12,
    )

    actual_electron_density = electron_density.set_axis(columns).to_frame(
        "value"
    )
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
    hydrogen_population_inputs: dict,
    regression_data: RegressionData,
) -> None:
    ion_population, electron_density, _ = solve_population(
        rate_matrix_solver, hydrogen_population_inputs, charge_conservation=True
    )
    shell = ion_population.columns[0]
    columns = pd.Index([0])
    actual_ion_population = ion_population[[shell]].set_axis(
        columns, axis="columns"
    )
    actual_ion_fraction = actual_ion_population / actual_ion_population.sum()
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
