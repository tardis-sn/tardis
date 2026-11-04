import numpy as np
import numpy.testing as npt
import pandas as pd
from astropy import units as u

from tardis.plasma.nonthermal.pynonthermal_solver import (
    NonthermalPlasmaInput,
    PynonthermalHeliumSolver,
)


class FakeSpencerFanoSolver:
    """Deterministic stand-in for the optional pynonthermal dependency."""

    def __init__(self, *args: object, **kwargs: object) -> None:
        del args, kwargs
        self.excitationlists = {(2, 1): {(0, 1): object()}}
        self.deposition_rate_density_ev = 0.0
        self.added_ions: list[tuple[int, int, float]] = []
        self.ionpopdict: dict[tuple[int, int], float] = {}

    def add_ionisation(self, Z: int, ion_stage: int, n_ion: float) -> None:
        self.added_ions.append((Z, ion_stage, n_ion))
        self.ionpopdict[(Z, ion_stage)] = n_ion

    def add_ion_ltepopexcitation(
        self,
        Z: int,
        ion_stage: int,
        n_ion: float,
        temperature: float,
    ) -> None:
        del Z, ion_stage, n_ion, temperature

    def solve(
        self, depositionratedensity_ev: float, override_n_e: float
    ) -> None:
        del override_n_e
        self.deposition_rate_density_ev = depositionratedensity_ev

    def analyse_ntspectrum(self) -> None:
        return None

    def get_excitation_ratecoeff(
        self, Z: int, ion_stage: int, transitionkey: tuple[int, int]
    ) -> float:
        del Z, ion_stage, transitionkey
        return 2.0

    def get_ionisation_ratecoeff(self, Z: int, ion_stage: int) -> float:
        return float(Z + ion_stage)

    def get_frac_heating(self) -> float:
        return 0.25


def make_input(deposition: u.Quantity) -> NonthermalPlasmaInput:
    """Return a two-shell helium plasma with aligned labels."""
    ion_index = pd.MultiIndex.from_tuples(
        [(2, 0), (2, 1), (2, 2)],
        names=["atomic_number", "ion_number"],
    )
    return NonthermalPlasmaInput(
        deposition,
        np.array([1.0e7, 2.0e7]) / u.cm**3,
        np.array([8000.0, 9000.0]) * u.K,
        pd.DataFrame(
            [[8.0e7, 7.0e7], [2.0e7, 3.0e7], [0.0, 0.0]],
            index=ion_index,
            columns=["inner", "outer"],
        ),
    )


def test_zero_deposition_returns_exact_zero_without_solver_construction() -> (
    None
):
    """A zero-input shell bypasses Spencer-Fano and has exact zero outputs."""

    def fail_factory(*args: object, **kwargs: object) -> object:
        del args, kwargs
        raise AssertionError("zero deposition constructed a solver")

    result = PynonthermalHeliumSolver(solver_factory=fail_factory).solve(
        make_input(np.zeros(2) * u.erg / u.s / u.cm**3)
    )

    npt.assert_array_equal(result.fixed_rates.bound_bound_rates, 0.0)
    npt.assert_array_equal(result.fixed_rates.ionization_rates, 0.0)
    npt.assert_array_equal(result.fixed_rates.heating_rate_density, 0.0)


def test_adapter_converts_units_stages_and_rate_normalizations() -> None:
    """Convert shell inputs and return fixed excitation, ionization, heating."""
    solvers: list[FakeSpencerFanoSolver] = []

    def factory(*args: object, **kwargs: object) -> FakeSpencerFanoSolver:
        solver = FakeSpencerFanoSolver(*args, **kwargs)
        solvers.append(solver)
        return solver

    deposition = np.array([1.0, 2.0]) * u.erg / u.s / u.cm**3
    result = PynonthermalHeliumSolver(solver_factory=factory).solve(
        make_input(deposition)
    )

    assert solvers[0].added_ions == [(2, 1, 8.0e7), (2, 2, 2.0e7)]
    excitation = result.fixed_rates.bound_bound_rates.loc[(2, 0, 0, 0, 0, 1)]
    npt.assert_allclose(
        excitation.to_numpy(),
        2.0 * deposition.to_value(u.eV / u.s / u.cm**3),
    )
    he_i_ionization = result.fixed_rates.ionization_rates.loc[(2, 0, 0, 1)]
    npt.assert_array_equal(he_i_ionization.to_numpy(), [3.0, 3.0])
    npt.assert_allclose(
        result.fixed_rates.heating_rate_density.to_numpy(), [0.25, 0.5]
    )


def test_superlevel_branching_preserves_total_excitation_rate() -> None:
    """Branching over detailed destinations conserves excitation rate."""
    destination_mapping = {1: ((4, 0.25), (5, 0.75))}
    result = PynonthermalHeliumSolver(
        destination_mapping=destination_mapping,
        solver_factory=FakeSpencerFanoSolver,
    ).solve(make_input(np.ones(2) * u.erg / u.s / u.cm**3))

    bound_rates = result.fixed_rates.bound_bound_rates
    npt.assert_allclose(
        bound_rates.loc[(2, 0, 0, 0, 0)].sum(axis=0),
        2.0 * (1.0 * u.erg).to_value(u.eV),
    )


def test_trace_ion_cutoff_preserves_required_helium_stages() -> None:
    """Discard trace competitors without removing helium rate providers.

    Claim: A composition cutoff removes ion stages below its fractional
    abundance while retaining positive He I and He II populations required by
    the adapter's requested rates.
    Regime: Two shells containing abundant helium and oxygen plus trace carbon.
    Verification: Inspect the independently recorded ions passed to the fake
    Spencer--Fano owner.
    """
    solvers: list[FakeSpencerFanoSolver] = []

    def factory(*args: object, **kwargs: object) -> FakeSpencerFanoSolver:
        solver = FakeSpencerFanoSolver(*args, **kwargs)
        solvers.append(solver)
        return solver

    plasma_input = make_input(np.ones(2) * u.erg / u.s / u.cm**3)
    plasma_input.ion_number_density.loc[(6, 0), :] = [1.0, 1.0]
    plasma_input.ion_number_density.loc[(8, 0), :] = [1.0e8, 1.0e8]

    PynonthermalHeliumSolver(
        solver_factory=factory,
        minimum_ion_number_fraction=1.0e-6,
    ).solve(plasma_input)

    for solver in solvers:
        added_stages = {
            (atomic_number, ion_stage)
            for atomic_number, ion_stage, _ in solver.added_ions
        }
        assert (2, 1) in added_stages
        assert (2, 2) in added_stages
        assert (8, 1) in added_stages
        assert (6, 1) not in added_stages
