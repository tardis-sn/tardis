from pathlib import Path

import pandas as pd
import pytest
from astropy import units as u

from tardis.io.atom_data import AtomData
from tardis.io.configuration.config_reader import Configuration
from tardis.model.base import SimulationState
from tardis.plasma.equilibrium.rate_matrix import AnalyticIonRateMatrix
from tardis.plasma.equilibrium.rates import (
    AnalyticPhotoionizationRateSolver,
    CollisionalIonizationRateSolver,
    RadiativeRatesSolver,
    ThermalCollisionalRateSolver,
)
from tardis.plasma.radiation_field.planck_rad_field import (
    DilutePlanckianRadiationField,
)

REFERENCE_ATOMIC_NUMBERS = (1, 2)
REFERENCE_SHELL_COUNT = 3
REFERENCE_ABUNDANCE = ((0.8, 0.8, 0.8), (0.2, 0.2, 0.2))
REFERENCE_DENSITIES = (1.0e-14, 1.1e-14, 1.2e-14)
REFERENCE_RADIATION_TEMPERATURES = (9500.0, 10000.0, 10500.0)
REFERENCE_DILUTION_FACTORS = (0.3, 0.5, 0.7)
REFERENCE_ELECTRON_TEMPERATURE_LINK = 0.9


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


@pytest.fixture
def new_chianti_atomic_dataset(tardis_regression_path):
    atomic_data_fname = (
        tardis_regression_path
        / "atom_data"
        / "kurucz_cd23_chianti_H_He_latest.h5"
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
def rate_solvers(
    radiative_rate_solver: RadiativeRatesSolver,
    collisional_rate_solver: ThermalCollisionalRateSolver,
) -> tuple[RadiativeRatesSolver, ThermalCollisionalRateSolver]:
    """Radiative and collisional rate solvers"""
    return radiative_rate_solver, collisional_rate_solver


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
def photoionization_rate_solver(mock_photoionization_cross_sections):
    return AnalyticPhotoionizationRateSolver(
        mock_photoionization_cross_sections
    )


@pytest.fixture
def collisional_ionization_rate_solver(mock_photoionization_cross_sections):
    return CollisionalIonizationRateSolver(mock_photoionization_cross_sections)


@pytest.fixture
def rate_matrix_solver(
    photoionization_rate_solver, collisional_ionization_rate_solver
):
    return AnalyticIonRateMatrix(
        photoionization_rate_solver, collisional_ionization_rate_solver
    )


@pytest.fixture
def mock_boltzmann_factor():
    index = pd.MultiIndex.from_tuples(
        [(1, 0, 0), (1, 0, 1)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    return pd.DataFrame([[2.0, 0.000011], [2.0, 0.003432]], index=index)


@pytest.fixture
def basic_thermodynamic_state(
    new_chianti_atomic_dataset: AtomData,
) -> dict:
    """Immutable, shared inputs for standard/IIP basic-state comparisons."""
    atomic_numbers = pd.Index(REFERENCE_ATOMIC_NUMBERS, name="atomic_number")
    columns = pd.Index(range(REFERENCE_SHELL_COUNT), name="shell")
    abundance = pd.DataFrame(
        REFERENCE_ABUNDANCE,
        index=atomic_numbers,
        columns=columns,
        dtype=float,
    )
    density = pd.Series(REFERENCE_DENSITIES, index=columns)
    t_rad = pd.Series(REFERENCE_RADIATION_TEMPERATURES, index=columns)
    dilution_factor = pd.Series(REFERENCE_DILUTION_FACTORS, index=columns)
    link_t_rad_t_electron = REFERENCE_ELECTRON_TEMPERATURE_LINK

    return {
        "atomic_data": new_chianti_atomic_dataset,
        "selected_atoms": atomic_numbers,
        "abundance": abundance,
        "density": density,
        "t_rad": t_rad,
        "dilution_factor": dilution_factor,
        "link_t_rad_t_electron": link_t_rad_t_electron,
        "radiation_field": DilutePlanckianRadiationField(
            t_rad.to_numpy() * u.K, dilution_factor.to_numpy()
        ),
    }
