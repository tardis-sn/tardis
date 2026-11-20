"""Case B must describe the same physical rows in level and ion equations."""

from types import MethodType

import numpy as np
import pandas as pd
import pytest
from astropy import units as u

from tardis import constants as const
from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.evaluator import PlasmaEquilibriumEvaluator
from tardis.plasma.equilibrium.inputs import (
    FixedRateState,
    ShellNumberDensity,
    SobolevInputs,
)
from tardis.plasma.equilibrium.ion_populations import IonPopulationSolver
from tardis.plasma.equilibrium.rate_matrix import (
    AnalyticIonRateMatrix,
    EstimatedIonRateMatrix,
    RateMatrix,
)
from tardis.plasma.equilibrium.rates import RadiativeRatesSolver
from tardis.plasma.equilibrium.rates.collisional_ionization_rates import (
    CollisionalIonizationRateSolver,
)
from tardis.plasma.equilibrium.rates.heating_cooling_rates import (
    BoundFreeThermalRates,
    CollisionalIonizationThermalRates,
)
from tardis.plasma.equilibrium.rates.photoionization_rates import (
    EstimatedPhotoionizationRateSolver,
)
from tardis.plasma.equilibrium.rates.util import (
    reindex_ionization_rate_dataframe,
    suppress_ground_recombination,
)
from tardis.plasma.equilibrium.solver import PlasmaEquilibriumSolver
from tardis.transport.montecarlo.estimators import init_estimators_continuum


@pytest.fixture
def helium_matrix_inputs() -> tuple:
    """Use the audit's exact atom; coefficient boundaries alone are prescribed."""
    index = pd.MultiIndex.from_tuples(
        [(2, 0, 0), (2, 0, 1), (2, 1, 0)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    photo = pd.DataFrame([2.0, 4.0, 5.0], index=index, columns=[18])
    spontaneous = pd.DataFrame([3.0, 5.0, 7.0], index=index, columns=[18])
    zero = photo * 0.0
    saha = zero + 1.0
    ions = pd.DataFrame(
        [4.0, 2.0, 1.0],
        columns=[18],
        index=pd.MultiIndex.from_tuples(
            [(2, 0), (2, 1), (2, 2)], names=["atomic_number", "ion_number"]
        ),
    )
    radiative = object.__new__(EstimatedPhotoionizationRateSolver)

    def coefficients(
        self: object, temperature: u.Quantity
    ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        return photo.copy(), zero.copy(), spontaneous.copy()

    radiative.solve_coefficients = MethodType(coefficients, radiative)
    collisional = CollisionalIonizationRateSolver(pd.DataFrame())
    collisional.solve_coefficients = lambda *args: zero.copy()
    fixed_index = pd.MultiIndex.from_tuples(
        [(2, 0, 0, 1)],
        names=[
            "atomic_number",
            "ion_number",
            "ion_number_source",
            "ion_number_destination",
        ],
    )
    owner = EstimatedIonRateMatrix(
        radiative,
        collisional,
        fixed_rates=FixedRateState(
            ionization_rates=pd.DataFrame(
                [1.0], index=fixed_index, columns=[18]
            )
        ),
    )
    distribution = ThermalElectronEnergyDistribution(
        0.0 * u.erg, np.array([8000.0]) * u.K, np.array([2.0]) / u.cm**3
    )
    return owner, distribution, photo, ions, saha


@pytest.mark.parametrize(
    "case_b, fractions, expected_ratio",
    [
        (False, [31 / 54, 23 / 54], 13 / 54),
        (True, [1 / 3, 2 / 3], 13 / 30),
    ],
)
def test_estimated_case_b_exact_ion_balance(
    helium_matrix_inputs: tuple,
    case_b: bool,
    fractions: list[float],
    expected_ratio: float,
) -> None:
    """The actual ion solve must reproduce the independently eliminated ratio."""
    owner, distribution, photo, ions, saha = helium_matrix_inputs
    owner.case_b_species = frozenset({(2, 0)} if case_b else ())
    levels = pd.DataFrame(
        np.r_[4 * np.array(fractions), 2.0],
        index=photo.index,
        columns=photo.columns,
    )
    owner.prepare(distribution, levels, ions, ions * 0 + 1, levels / 4, saha)
    matrix = owner.solve_prepared(np.array([2.0]), photo.columns).at[2, 18]
    population = np.linalg.solve(matrix, np.array([0.0, 1.0, 0.0]))
    np.testing.assert_allclose(
        population[1] / population[0], expected_ratio, rtol=1e-13
    )


@pytest.mark.parametrize(
    "owner_type", [AnalyticIonRateMatrix, EstimatedIonRateMatrix]
)
@pytest.mark.parametrize("electron_density", [1.0, 3.0, 7.0])
def test_case_b_matrix_density_powers_and_labels(
    owner_type: type, electron_density: float
) -> None:
    """Mask only He I ground: recombination is ne*(spont+stim)+ne²*threebody.

    Deliberately non-positional ground rows and a He II ground distinguish
    bound destinations from continuum source ions. Hydrogen remains unchanged.
    Only coefficient boundaries are prescribed; all matrix assembly is real.
    """
    index = pd.MultiIndex.from_tuples(
        [(2, 0, 7), (1, 0, 0), (2, 1, 0), (2, 0, 0)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    spontaneous = pd.DataFrame([3.0, 5.0, 7.0, 11.0], index=index, columns=[18])
    stimulated = spontaneous * 2
    threebody = spontaneous * 5
    photo = spontaneous * 7
    coll_ion = spontaneous * 11
    originals = [
        frame.copy(deep=True)
        for frame in (spontaneous, stimulated, threebody, photo, coll_ion)
    ]

    class Radiative:
        def solve(self, *args: object) -> tuple[pd.DataFrame, pd.DataFrame]:
            density = (
                electron_density if owner_type is AnalyticIonRateMatrix else 1.0
            )
            return reindex_ionization_rate_dataframe(
                photo.copy()
            ), reindex_ionization_rate_dataframe(
                (spontaneous + stimulated) * density, recombination=True
            )

    class Collisional:
        def solve(
            self, *args: object, **kwargs: object
        ) -> tuple[pd.DataFrame, pd.DataFrame]:
            density = (
                electron_density if owner_type is AnalyticIonRateMatrix else 1.0
            )
            return reindex_ionization_rate_dataframe(
                coll_ion * density
            ), reindex_ionization_rate_dataframe(
                threebody * density**2, recombination=True
            )

    owner = owner_type(Radiative(), Collisional())
    distribution = ThermalElectronEnergyDistribution(
        0.0 * u.erg,
        np.array([8000.0]) * u.K,
        np.array([electron_density]) / u.cm**3,
    )

    def matrix_solve() -> pd.DataFrame:
        return owner.solve(
            None,
            distribution,
            photo,
            photo,
            photo,
            photo,
            photo,
            photo,
            spontaneous,
        )

    owner.case_b_species = ()
    case_a = matrix_solve()
    owner.case_b_species = ((2, 0),)
    case_b = matrix_solve()
    np.testing.assert_allclose(
        case_a.at[2, 18][0, 1],
        (3 + 11) * (3 * electron_density + 5 * electron_density**2),
        rtol=1e-14,
    )
    np.testing.assert_allclose(
        case_b.at[2, 18][0, 1],
        3 * (3 * electron_density + 5 * electron_density**2),
        rtol=1e-14,
    )
    np.testing.assert_allclose(
        case_a.at[2, 18][0, 0], case_b.at[2, 18][0, 0], rtol=0, atol=0
    )
    np.testing.assert_allclose(
        case_a.at[2, 18][2, 2], case_b.at[2, 18][2, 2], rtol=0, atol=0
    )
    np.testing.assert_array_equal(case_a.at[1, 18], case_b.at[1, 18])
    owner.case_b_species = ()
    restored = matrix_solve()
    np.testing.assert_array_equal(restored.at[2, 18], case_a.at[2, 18])
    for original, frame in zip(
        originals,
        (spontaneous, stimulated, threebody, photo, coll_ion),
        strict=True,
    ):
        pd.testing.assert_frame_equal(frame, original)


def test_estimated_policy_invalidates_prepared_cache(
    helium_matrix_inputs: tuple,
) -> None:
    """Same policy preserves cache; B→A→B discards every stale polynomial."""
    owner, distribution, photo, ions, saha = helium_matrix_inputs
    owner.prepare(distribution, photo, ions, ions * 0 + 1, photo, saha)
    prepared = owner.prepared_matrix_coefficients
    owner.case_b_species = ()
    assert owner.prepared_matrix_coefficients is prepared
    for policy in [((2, 0),), (), ((2, 0),)]:
        owner.case_b_species = policy
        assert owner.prepared_matrix_coefficients is None
        assert owner.prepared_columns is None
        owner.prepare(distribution, photo, ions, ions * 0 + 1, photo, saha)


def test_ground_mask_is_copying_and_idempotent() -> None:
    """Mask actual level zero, not first row or continuum's ion label."""
    index = pd.MultiIndex.from_tuples(
        [(2, 0, 5), (2, 0, 0), (2, 1, 0), (1, 0, 0)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    frame = pd.DataFrame([1.0, 2.0, 3.0, 4.0], index=index, columns=[53])
    original = frame.copy(deep=True)
    masked = suppress_ground_recombination(frame, ((2, 0),))
    np.testing.assert_array_equal(masked[53], [1.0, 0.0, 3.0, 4.0])
    pd.testing.assert_frame_equal(
        masked, suppress_ground_recombination(masked, ((2, 0),))
    )
    pd.testing.assert_frame_equal(frame, original)


def test_explicit_case_b_thermal_terms() -> None:
    """Remove He I ground recombination thermal terms, never ionization losses."""
    index = pd.MultiIndex.from_tuples(
        [(1, 0, 0), (2, 0, 0), (2, 0, 1), (2, 1, 0)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    sections = pd.DataFrame(
        {"nu": np.tile([1.0e15, 2.0e15], 4), "x_sect": np.ones(8) * 1.0e-18},
        index=index.repeat(2),
    )
    ions = pd.DataFrame(
        [5.0, 7.0, 11.0, 13.0, 17.0],
        index=pd.MultiIndex.from_tuples(
            [(1, 0), (1, 1), (2, 0), (2, 1), (2, 2)],
            names=["atomic_number", "ion_number"],
        ),
        columns=[18],
    )
    levels = pd.DataFrame(1.0, index=index, columns=[18])
    distribution = ThermalElectronEnergyDistribution(
        0.0 * u.erg, np.array([8000.0]) * u.K, np.array([2.0]) / u.cm**3
    )
    stimulated = levels * 3.0
    free_bound = BoundFreeThermalRates(sections)
    legacy = free_bound.calculate_cooling_rate(ions, distribution, levels)
    legacy_stim = free_bound.calculate_cooling_rate(
        ions, distribution, levels, stimulated
    )
    case_b = free_bound.calculate_cooling_rate(
        ions, distribution, levels, stimulated, ((2, 0),)
    )
    assert legacy.loc[(1, 0, 0), 18] == 0.0
    assert (
        legacy_stim.loc[(1, 0, 0), 18] > 0.0
    )  # unchanged legacy H stimulated channel
    expected = legacy_stim.copy()
    expected.loc[(2, 0, 0)] = 0.0
    pd.testing.assert_frame_equal(case_b, expected)
    heating_estimator = levels * 5.0
    heating_a, _ = free_bound.solve(
        levels,
        ions,
        distribution,
        levels,
        bound_free_heating_estimator=heating_estimator,
        stimulated_recombination_estimator=stimulated,
    )
    heating_b, _ = free_bound.solve(
        levels,
        ions,
        distribution,
        levels,
        bound_free_heating_estimator=heating_estimator,
        stimulated_recombination_estimator=stimulated,
        case_b_species=((2, 0),),
    )
    pd.testing.assert_series_equal(heating_a, heating_b)
    collisional = CollisionalIonizationThermalRates(sections)
    heating_a, cooling_a = collisional.solve(
        distribution.number_density, ions, levels, levels, levels
    )
    heating_b, cooling_b = collisional.solve(
        distribution.number_density, ions, levels, levels, levels, ((2, 0),)
    )
    np.testing.assert_allclose(
        (heating_a - heating_b).iloc[0],
        2.0**2 * 13.0 * 1.0e15 * const.h.cgs.value,
        rtol=1e-14,
    )
    pd.testing.assert_series_equal(cooling_a, cooling_b)
    pd.testing.assert_frame_equal(stimulated, levels * 3.0)


@pytest.mark.parametrize("policy", [(), ((2, 0),)])
def test_evaluator_case_b_terminal_charge_and_level_closure(
    helium_matrix_inputs: tuple, policy: tuple
) -> None:
    """Actual evaluator/ion/outer solves close one synthetic helium plasma.

    Prescribe positive continuum coefficients and negligible line optical depth;
    the density-temperature Saha conversion remains real. Check conservation and
    terminal kinetic rows independently of optimizer success flags. No helium
    policy is duplicated on the ion owner, which is the reported wiring failure.
    """
    owner, _distribution, photo, ions, _ = helium_matrix_inputs
    bound_index = photo.index
    level_index = bound_index.append(
        pd.MultiIndex.from_tuples([(2, 2, 0)], names=bound_index.names)
    )
    levels = pd.DataFrame(
        {"energy": [0.0, 1.0e-12, 0.0, 0.0], "g": [1.0, 3.0, 2.0, 1.0]},
        index=level_index,
    )
    lines = pd.MultiIndex.from_tuples(
        [(2, 0, 0, 1)],
        names=[
            "atomic_number",
            "ion_number",
            "level_number_lower",
            "level_number_upper",
        ],
    )

    class ElectronRates:
        all_collisional_strengths_index = pd.MultiIndex.from_tuples(
            [],
            names=[
                "atomic_number",
                "ion_number",
                "ion_number_source",
                "ion_number_destination",
                "level_number_source",
                "level_number_destination",
            ],
        )

        def solve(self, temperature: u.Quantity) -> pd.DataFrame:
            return pd.DataFrame(
                index=self.all_collisional_strengths_index,
                columns=[18],
                dtype=float,
            )

    bound_owner = RateMatrix(
        RadiativeRatesSolver(
            pd.DataFrame(
                {"A_ul": [2.0], "B_ul": [0.0], "B_lu": [1.0], "nu": [1.0e15]},
                index=lines,
            )
        ),
        ElectronRates(),
        levels,
    )
    evaluator = PlasmaEquilibriumEvaluator(
        pd.DataFrame(
            {"nu": [1.0e15, 1.1e15, 1.2e15], "x_sect": [1.0e-18] * 3},
            index=bound_index,
        ),
        pd.Series(np.arange(3), index=bound_index),
        init_estimators_continuum((3, 1), 1),
        1.0 * u.s,
        1.0 * u.cm**3,
        levels,
        pd.Series(
            [2.0e-11, 4.0e-11],
            index=pd.MultiIndex.from_tuples(
                [(2, 1), (2, 2)], names=["atomic_number", "ion_number"]
            ),
        ),
        bound_owner,
        pd.DataFrame([1.0], index=lines, columns=[18]),
        (
            ShellNumberDensity(
                7.0, np.array([2.0, 2.0, 2.0, 1.0]), np.array([0, 1])
            ),
        ),
        (
            SobolevInputs(
                np.array([0]),
                np.array([1]),
                np.array([1.0]),
                np.array([3.0]),
                np.array([False]),
                np.array([False]),
                np.array([0.0]),
                np.array([0]),
                lines,
            ),
        ),
        level_index,
        (2, 0),
        pd.DataFrame([7.0], index=[2], columns=[18]),
        np.array([14.0]),
        ion_population_solver=IonPopulationSolver(owner),
        ion_population_arguments={
            "radiation_field": None,
            "elemental_number_density": pd.DataFrame(
                [7.0], index=[2], columns=[18]
            ),
            "estimated_ion_population": ions,
        },
        fixed_rates=owner.fixed_rates,
        case_b_species=policy,
    )
    saha = evaluator.calculate_continuum_coefficients(np.array([8000.0]))[1]
    original_coefficients = (
        owner.radiative_ionization_rate_solver.solve_coefficients
    )

    def coefficients(
        self: object, temperature: u.Quantity
    ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        photo_coeff, stim_coeff, spontaneous_coeff = original_coefficients(
            temperature
        )
        return photo_coeff, stim_coeff, spontaneous_coeff / saha

    owner.radiative_ionization_rate_solver.solve_coefficients = MethodType(
        coefficients, owner.radiative_ionization_rate_solver
    )
    evaluator.estimated_photoionization_rate_solver = (
        owner.radiative_ionization_rate_solver
    )
    evaluator.collisional_ionization_rate_solver = (
        owner.collisional_ionization_rate_solver
    )
    seed = pd.DataFrame([0.5, 0.5], index=bound_index[:2], columns=[18])
    solution = PlasmaEquilibriumSolver(fixed_temperature_ratio=1.0).solve(
        evaluator, [2.0], [8000.0], seed, [14.0], [8000.0]
    )
    result = solution.evaluation
    assert owner.case_b_species == frozenset(policy)
    np.testing.assert_allclose(result.ion_population.sum(), [7.0], rtol=1.0e-12)
    charge = result.ion_population.multiply(
        result.ion_population.index.get_level_values("ion_number"), axis=0
    ).sum()
    np.testing.assert_allclose(charge, solution.electron_density, rtol=1.0e-9)
    totals = result.absolute_level_population.groupby(
        level=["atomic_number", "ion_number"]
    ).sum()
    pd.testing.assert_frame_equal(
        totals, result.ion_population.loc[totals.index], rtol=1.0e-12
    )
    assert np.max(np.abs(result.relative_level_residual.to_numpy())) < 1.0e-8
    neutral = result.absolute_level_population.loc[(2, 0), 18].to_numpy()
    ionized = result.ion_population.loc[(2, 1), 18]
    density = solution.electron_density[0]
    recombination = density * np.array([0.0 if policy else 3.0, 5.0])
    bound_flows = np.array(
        [-neutral[0] + 2 * neutral[1], neutral[0] - 2 * neutral[1]]
    )
    losses = np.array([3.0, 5.0]) * neutral
    sources = recombination * ionized
    np.testing.assert_allclose(
        (bound_flows - losses + sources)
        / (np.abs(bound_flows) + losses + sources),
        0.0,
        atol=1.0e-8,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        result.ion_population.loc[(2, 2), 18] / ionized,
        5.0 / (7.0 * density),
        rtol=1.0e-8,
    )
