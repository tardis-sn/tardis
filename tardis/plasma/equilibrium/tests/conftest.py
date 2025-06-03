from pathlib import Path

import pandas as pd
import pytest

from tardis.io.atom_data import AtomData
from tardis.io.configuration.config_reader import Configuration
from tardis.model.base import SimulationState
from tardis.plasma.equilibrium.rate_matrix import IonRateMatrix
from tardis.plasma.equilibrium.rates import (
    AnalyticPhotoionizationRateSolver,
    CollisionalIonizationRateSolver,
    RadiativeRatesSolver,
    ThermalCollisionalRateSolver,
)


@pytest.fixture
def mock_photoionization_cross_sections():
    """Fixture for mock photoionization cross-sections."""
    data = {
        "nu": [1e15, 2e15],
        "x_sect": [1e-18, 2e-18],
    }
    index = pd.MultiIndex.from_tuples(
        [(1, 0, 0), (1, 0, 1)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    return pd.DataFrame(data, index=index)


@pytest.fixture
def new_chianti_atomic_dataset_si(tardis_regression_path):
    atomic_data_fname = (
        tardis_regression_path / "atom_data" / "kurucz_cd23_chianti_Si.h5"
    )
    return AtomData.from_hdf(atomic_data_fname)


@pytest.fixture(scope="session")
def hydrogen_atomic_data_fname(tardis_regression_path):
    """
    File name for the atomic data file used in NTLE ionization solver tests.
    """
    atomic_data_fname = (
        tardis_regression_path / "atom_data" / "nlte_atom_data" / "cmfgen_H.h5"
    )

    atom_data_missing_str = (
        f"{atomic_data_fname} atomic datafiles does not seem to exist"
    )

    if not atomic_data_fname.exists():
        pytest.exit(atom_data_missing_str)

    return atomic_data_fname


@pytest.fixture(scope="session")
def hydrogen_atomic_dataset(hydrogen_atomic_data_fname):
    """
    Atomic dataset used for NLTE ionization solver tests.
    """
    h_atomic_data = AtomData.from_hdf(hydrogen_atomic_data_fname)
    return h_atomic_data


@pytest.fixture(params=[(14, 1, slice(None), slice(None))])
def radiative_transitions(new_chianti_atomic_dataset_si, request):
    return new_chianti_atomic_dataset_si.lines.loc[request.param, :]


@pytest.fixture
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


@pytest.fixture
def rate_solver_list(radiative_rate_solver, collisional_rate_solver):
    return [
        (radiative_rate_solver, "radiative"),
        (collisional_rate_solver, "electron"),
    ]


@pytest.fixture
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


@pytest.fixture
def photoionization_rate_solver(hydrogen_atomic_dataset):
    return AnalyticPhotoionizationRateSolver(
        hydrogen_atomic_dataset.photoionization_data
    )


@pytest.fixture
def collisional_ionization_rate_solver(hydrogen_atomic_dataset):
    return CollisionalIonizationRateSolver(
        hydrogen_atomic_dataset.photoionization_data
    )


@pytest.fixture
def rate_matrix_solver(
    photoionization_rate_solver, collisional_ionization_rate_solver
):
    return IonRateMatrix(
        photoionization_rate_solver, collisional_ionization_rate_solver
    )
