from pathlib import Path

import pytest

from tardis.io.atom_data import AtomData
from tardis.io.configuration.config_reader import Configuration
from tardis.model.base import SimulationState
from tardis.plasma.equilibrium.rates import (
    RadiativeRatesSolver,
    ThermalCollisionalRateSolver,
)


@pytest.fixture()
def new_chianti_atomic_dataset_si(tardis_regression_path):
    atomic_data_fname = (
        tardis_regression_path / "atom_data" / "kurucz_cd23_chianti_Si.h5"
    )
    return AtomData.from_hdf(atomic_data_fname)


@pytest.fixture(params=[(14, 1, slice(None), slice(None))])
def radiative_transitions(new_chianti_atomic_dataset_si, request):
    return new_chianti_atomic_dataset_si.lines.loc[request.param, :]


@pytest.fixture()
def radiative_rate_solver(radiative_transitions):
    return RadiativeRatesSolver(radiative_transitions)


@pytest.fixture(params=[(14, 1, slice(None), slice(None))])
def collisional_rate_solver(
    new_chianti_atomic_dataset_si, radiative_transitions, request
):
    col_strength_temperatures = (
        new_chianti_atomic_dataset_si.collision_data_temperatures
    )
    col_strengths = new_chianti_atomic_dataset_si.collision_data.loc[
        request.param, :
    ]
    return ThermalCollisionalRateSolver(
        new_chianti_atomic_dataset_si.levels,
        radiative_transitions,
        col_strength_temperatures,
        col_strengths,
        "chianti",
    )


@pytest.fixture()
def rate_solver_list(radiative_rate_solver, collisional_rate_solver):
    return [
        (radiative_rate_solver, "radiative"),
        (collisional_rate_solver, "electron"),
    ]


@pytest.fixture()
def collisional_simulation_state(new_chianti_atomic_dataset_si):
    config = Configuration.from_yaml(
        Path("tardis")
        / "plasma"
        / "tests"
        / "data"
        / "plasma_base_test_config.yml"
    )
    return SimulationState.from_config(
        config, atom_data=new_chianti_atomic_dataset_si
    )
