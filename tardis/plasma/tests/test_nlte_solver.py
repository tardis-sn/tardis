import os

import pytest
import numpy as np
import pandas as pd
from copy import deepcopy
from numpy.testing import assert_allclose, assert_almost_equal
from tardis.io.config_reader import Configuration
from tardis.model.base import Radial1DModel
from tardis.plasma.properties import NLTERateEquationSolver
from tardis.io.atom_data.base import AtomData
from tardis.plasma.standard_plasmas import assemble_plasma


@pytest.fixture
def simple_index_nlte_ion():
    """Simple fixture for nlte_ion treatment for H I and He II.

    Returns
    -------
    MultiIndex
        MultiIndex for HI and HeII.
    """
    return pd.MultiIndex.from_tuples(
        [(1, 0), (2, 1)], names=("atomic_number", "ion_number")
    )


@pytest.fixture
def simple_index_lte_ion():
    """Simple fixture for lte_ion treatment for H II and He I and He III.

    Returns
    -------
    MultiIndex
        MultiIndex for H II and He I and He III.
    """
    return pd.MultiIndex.from_tuples(
        [(1, 1), (2, 1), (2, 2)], names=("atomic_number", "ion_number")
    )


@pytest.fixture
def simple_rate_matrix_index():
    """Simple rate_matrix_index for NTLE ionization treatment of H I and He II.

    Returns
    -------
    MultiIndex
        (atomic_number, ion_number, treatment)
    """
    return pd.MultiIndex.from_tuples(
        [
            (1, 0, "nlte_ion"),
            (1, 1, "lte_ion"),
            (2, 0, "lte_ion"),
            (2, 1, "nlte_ion"),
            (2, 2, "lte_ion"),
            ("n_e", "n_e", "n_e"),
        ],
        names=("atomic_number", "ion_number", "level_number"),
    )


@pytest.fixture
def simple_total_photo_ion_coefficients(simple_index_nlte_ion):
    """Simple coefficients for photoionization of H I and He II.

    Returns
    -------
    DataFrame
        Photoionization coefficients for H I and He II.
    """
    simple_photo_ion_coefficients = [0.03464792, 0.68099508]
    return pd.DataFrame(
        simple_photo_ion_coefficients, index=simple_index_nlte_ion
    )


@pytest.fixture
def simple_total_rad_recomb_coefficients(simple_index_nlte_ion):
    """Simple coefficients for radiative recombination of H I and He II.

    Returns
    -------
    DataFrame
        Radiative recombination coefficients for H I and He II.
    """
    simple_rad_recomb_coefficients = [0.43303813, 0.66140309]
    return pd.DataFrame(
        simple_rad_recomb_coefficients, index=simple_index_nlte_ion
    )


@pytest.fixture
def simple_total_col_ion_coefficients(simple_index_nlte_ion):
    """Simple coefficients for collisional ionization of H I and He II.

    Returns
    -------
    DataFrame
        Collisional ionization coefficients for H I and He II.
    """
    simple_col_ion_coefficients = [0.19351674, 0.69214007]
    return pd.DataFrame(
        simple_col_ion_coefficients, index=simple_index_nlte_ion
    )


@pytest.fixture
def simple_total_col_recomb_coefficients(simple_index_nlte_ion):
    """Simple coefficients for collisional recombination of H I and He II.

    Returns
    -------
    DataFrame
        Collisional recombination coefficients for H I and He II.
    """
    simple_col_recomb_coefficients = [0.06402515, 0.29785023]
    return pd.DataFrame(
        simple_col_recomb_coefficients, index=simple_index_nlte_ion
    )


@pytest.fixture
def simple_phi(simple_index_lte_ion):
    """Simple Saha factors for H II, He I and He III."""
    simple_phi = [0.18936306, 0.15726292, 0.79851244]
    return pd.DataFrame(simple_phi, index=simple_index_lte_ion)


@pytest.fixture
def simple_electron_density():
    """Simple electron density."""
    return 0.2219604493076


def test_rate_matrix(
    simple_phi,
    simple_electron_density,
    simple_rate_matrix_index,
    simple_total_photo_ion_coefficients,
    simple_total_rad_recomb_coefficients,
    simple_total_col_ion_coefficients,
    simple_total_col_recomb_coefficients,
):
    """
    Using a simple case of nlte_ion for HI and HeII, checks if the calculate_rate_matrix generates the correct data.
    """
    atomic_numbers = [1, 2]
    actual_rate_matrix = NLTERateEquationSolver.calculate_rate_matrix(
        atomic_numbers,
        simple_phi,
        simple_electron_density,
        simple_rate_matrix_index,
        simple_total_photo_ion_coefficients,
        simple_total_rad_recomb_coefficients,
        simple_total_col_ion_coefficients,
        simple_total_col_recomb_coefficients,
    )
    desired_rate_matrix = [
        [-0.077601, 0.099272, 0.000000, 0.000000, 0.000000, 0.0],
        [1.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.0],
        [0.000000, 0.000000, -0.157263, 0.221960, 0.000000, 0.0],
        [0.000000, 0.000000, 0.000000, -0.834623, 0.161479, 0.0],
        [0.000000, 0.000000, 1.000000, 1.000000, 1.000000, 0.0],
        [0.000000, 1.000000, 0.000000, 1.000000, 2.000000, -1.0],
    ]

    assert_almost_equal(
        desired_rate_matrix, np.array(actual_rate_matrix), decimal=6
    )


def test_jacobian_matrix(
    simple_phi,
    simple_electron_density,
    simple_rate_matrix_index,
    simple_total_photo_ion_coefficients,
    simple_total_rad_recomb_coefficients,
    simple_total_col_ion_coefficients,
    simple_total_col_recomb_coefficients,
):
    """
    Using a simple case of nlte_ion for HI and HeII,
    checks if the jacobian_matrix generates the correct data.
    """
    atomic_numbers = [1, 2]

    initial_guess = [
        0.7192433675307516,
        0.8101666197902874,
        0.7171853313284426,
        0.040220760173800496,
        0.2878574499274399,
        simple_electron_density,
    ]
    simple_rate_matrix = NLTERateEquationSolver.calculate_rate_matrix(
        atomic_numbers,
        simple_phi,
        simple_electron_density,
        simple_rate_matrix_index,
        simple_total_photo_ion_coefficients,
        simple_total_rad_recomb_coefficients,
        simple_total_col_ion_coefficients,
        simple_total_col_recomb_coefficients,
    )

    actual_jacobian_matrix = NLTERateEquationSolver.jacobian_matrix(
        atomic_numbers,
        initial_guess,
        simple_rate_matrix,
        simple_rate_matrix_index,
        simple_total_rad_recomb_coefficients,
        simple_total_col_ion_coefficients,
        simple_total_col_recomb_coefficients,
    )

    desired_jacobian_matrix = [
        [-0.07760098, 0.09927163, 0.0, 0.0, 0.0, 0.23467404],
        [1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, -0.15726292, 0.22196045, 0.0, 0.04022076],
        [0.0, 0.0, 0.0, -0.8346228, 0.16147935, 0.20061248],
        [0.0, 0.0, 1.0, 1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0, 1.0, 2.0, -1.0],
    ]
    assert_almost_equal(actual_jacobian_matrix, desired_jacobian_matrix)


@pytest.fixture(scope="session")
def nlte_atomic_data_fname(tardis_ref_path):
    atomic_data_fname = os.path.join(
        tardis_ref_path, "nlte_atom_data", "TestNLTE_He_Ti.h5"
    )

    atom_data_missing_str = (
        f"{atomic_data_fname} atomic datafiles " f"does not seem to exist"
    )

    if not os.path.exists(atomic_data_fname):
        pytest.exit(atom_data_missing_str)

    return atomic_data_fname


@pytest.fixture(scope="session")
def nlte_atomic_dataset(nlte_atomic_data_fname):
    nlte_atomic_data = AtomData.from_hdf(nlte_atomic_data_fname)

    # if atomic_data.md5 != DEFAULT_ATOM_DATA_UUID:
    #     pytest.skip(
    #         f'Need default Kurucz atomic dataset (md5="{DEFAULT_ATOM_DATA_UUID}")'
    #     )
    # else:
    return nlte_atomic_data


@pytest.fixture
def nlte_atom_data(nlte_atomic_dataset):
    atomic_data = deepcopy(nlte_atomic_dataset)
    return atomic_data


data_path = os.path.join("tardis", "io", "tests", "data")


@pytest.fixture
def tardis_model_config_nlte():
    filename = "tardis_configv1_nlte.yml"
    return Configuration.from_yaml(os.path.join(data_path, filename))


@pytest.fixture
def nlte_raw_model(tardis_model_config_nlte):
    return Radial1DModel.from_config(tardis_model_config_nlte)


@pytest.fixture
def nlte_raw_plasma(tardis_model_config_nlte, nlte_raw_model, nlte_atom_data):
    new_w = np.ones_like(nlte_raw_model.dilution_factor)
    nlte_raw_model.dilution_factor = new_w
    plasma = assemble_plasma(
        tardis_model_config_nlte, nlte_raw_model, nlte_atom_data
    )
    return plasma


def test_critical_case(nlte_raw_plasma):
    ion_number_density_nlte = nlte_raw_plasma.ion_number_density_nlte.values
    ion_number_density_nlte[ion_number_density_nlte < 1e-10] = 0.0

    ion_number_density = nlte_raw_plasma.ion_number_density.values
    ion_number_density[ion_number_density < 1e-10] = 0.0
    assert_allclose(
        ion_number_density,
        ion_number_density_nlte,
        rtol=1e-2 * ion_number_density.max(),
    )
