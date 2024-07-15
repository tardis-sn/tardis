import numpy as np
import numpy.testing as npt
import pandas as pd
import pytest
import numpy.testing as npt


from tardis.plasma.properties import (
    NLTEPopulationSolverLU,
    NLTEPopulationSolverRoot,
)
from tardis.plasma.properties.ion_population import IonNumberDensity
from tardis.plasma.properties.nlte_rate_equation_solver import (
    calculate_jacobian_matrix,
    calculate_rate_matrix,
)
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
    regression_data,
):
    """
    Using a simple case of nlte_ion for HI and HeII, checks if the calculate_rate_matrix generates the correct data.
    """
    atomic_numbers = [1, 2]
    actual_rate_matrix = calculate_rate_matrix(
        atomic_numbers,
        simple_phi,
        simple_electron_density,
        simple_rate_matrix_index,
        simple_total_photo_ion_coefficients,
        simple_total_rad_recomb_coefficients,
        simple_total_col_ion_coefficients,
        simple_total_col_recomb_coefficients,
    )
    # TODO: decimal=6
    # allow for assert_almost_equal
    expected_rate_matrix = regression_data.sync_ndarray(actual_rate_matrix)
    npt.assert_allclose(actual_rate_matrix, expected_rate_matrix, rtol=1e-6)


def test_jacobian_matrix(
    simple_phi,
    simple_electron_density,
    simple_rate_matrix_index,
    simple_total_photo_ion_coefficients,
    simple_total_rad_recomb_coefficients,
    simple_total_col_ion_coefficients,
    simple_total_col_recomb_coefficients,
    regression_data,
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
    simple_rate_matrix = calculate_rate_matrix(
        atomic_numbers,
        simple_phi,
        simple_electron_density,
        simple_rate_matrix_index,
        simple_total_photo_ion_coefficients,
        simple_total_rad_recomb_coefficients,
        simple_total_col_ion_coefficients,
        simple_total_col_recomb_coefficients,
    )

    actual_jacobian_matrix = calculate_jacobian_matrix(
        atomic_numbers,
        initial_guess,
        simple_rate_matrix,
        simple_rate_matrix_index,
        simple_total_rad_recomb_coefficients,
        simple_total_col_ion_coefficients,
        simple_total_col_recomb_coefficients,
    )

    # TODO: allow for assert_almost_equal
    expected_jacobian_matrix = regression_data.sync_ndarray(
        actual_jacobian_matrix
    )
    npt.assert_allclose(actual_jacobian_matrix, expected_jacobian_matrix)


@pytest.fixture
def nlte_raw_plasma_dilution_factor_1_root(
    tardis_model_config_nlte_root, nlte_raw_model_root, nlte_atom_data
):
    """
    Plasma assembled with dilution factors set to 1.0.
    """
    new_dilution_factor = np.ones_like(nlte_raw_model_root.dilution_factor)
    nlte_raw_model_root.dilution_factor = new_dilution_factor
    plasma = assemble_plasma(
        tardis_model_config_nlte_root, nlte_raw_model_root, nlte_atom_data
    )
    return plasma


@pytest.fixture
def nlte_raw_plasma_dilution_factor_1_lu(
    tardis_model_config_nlte_lu, nlte_raw_model_lu, nlte_atom_data
):
    """
    Plasma assembled with dilution factors set to 1.0.
    """
    new_dilution_factor = np.ones_like(nlte_raw_model_lu.dilution_factor)
    nlte_raw_model_lu.dilution_factor = new_dilution_factor
    plasma = assemble_plasma(
        tardis_model_config_nlte_lu, nlte_raw_model_lu, nlte_atom_data
    )
    return plasma


@pytest.fixture
def nlte_raw_plasma_dilution_factor_0_root(
    tardis_model_config_nlte_root, nlte_raw_model_root, nlte_atom_data
):
    """
    Plasma assembled with dilution factors set to 0.0.
    """
    new_dilution_factor = np.zeros_like(nlte_raw_model_root.dilution_factor)
    nlte_raw_model_root.dilution_factor = new_dilution_factor
    plasma = assemble_plasma(
        tardis_model_config_nlte_root, nlte_raw_model_root, nlte_atom_data
    )
    return plasma


@pytest.fixture
def nlte_raw_plasma_dilution_factor_0_lu(
    tardis_model_config_nlte_lu, nlte_raw_model_lu, nlte_atom_data
):
    """
    Plasma assembled with dilution factors set to 0.0.
    """
    new_dilution_factor = np.zeros_like(nlte_raw_model_lu.dilution_factor)
    nlte_raw_model_lu.dilution_factor = new_dilution_factor
    plasma = assemble_plasma(
        tardis_model_config_nlte_lu, nlte_raw_model_lu, nlte_atom_data
    )
    return plasma


def test_critical_case_dilution_factor_1_root(
    nlte_raw_plasma_dilution_factor_1_root,
):
    """Check that the LTE and NLTE solution agree for dilution_factor=1.0."""
    ion_number_density_nlte = (
        nlte_raw_plasma_dilution_factor_1_root.ion_number_density.values
    )
    # ion_number_density_nlte[np.abs(ion_number_density_nlte) < 1e-10] = 0.0

    ind = IonNumberDensity(nlte_raw_plasma_dilution_factor_1_root)
    ion_number_density_lte = ind.calculate(
        nlte_raw_plasma_dilution_factor_1_root.thermal_phi_lte,
        nlte_raw_plasma_dilution_factor_1_root.partition_function,
        nlte_raw_plasma_dilution_factor_1_root.number_density,
    )[0]

    npt.assert_allclose(
        ion_number_density_lte,
        ion_number_density_nlte,
        atol=1e-10,  # seems fair for the test
        rtol=1e-2,
    )


def test_critical_case_dilution_factor_1_lu(
    nlte_raw_plasma_dilution_factor_1_lu,
):
    """Check that the LTE and NLTE solution agree for dilution_factor=1.0."""
    ion_number_density_nlte = (
        nlte_raw_plasma_dilution_factor_1_lu.ion_number_density.values
    )
    # ion_number_density_nlte[np.abs(ion_number_density_nlte) < 1e-10] = 0.0

    ind = IonNumberDensity(nlte_raw_plasma_dilution_factor_1_lu)
    ion_number_density_lte = ind.calculate(
        nlte_raw_plasma_dilution_factor_1_lu.thermal_phi_lte,
        nlte_raw_plasma_dilution_factor_1_lu.partition_function,
        nlte_raw_plasma_dilution_factor_1_lu.number_density,
    )[0]

    ion_number_density_lte = ion_number_density_lte.values
    ion_number_density_lte[
        ion_number_density_lte < 1e-10
    ] = 0.0  # getting rid of small numbers.
    npt.assert_allclose(
        ion_number_density_lte,
        ion_number_density_nlte,
        atol=1e-10,  # seems fair for the test
        rtol=1e-2,
    )


def test_critical_case_dilution_factor_0_root(
    nlte_raw_plasma_dilution_factor_0_root,
):
    """Check that the LTE and NLTE solution agree for dilution_factor=0.0."""
    nlte_solver = NLTEPopulationSolverRoot(
        nlte_raw_plasma_dilution_factor_0_root
    )
    ion_number_density_nlte = nlte_solver.calculate(
        nlte_raw_plasma_dilution_factor_0_root.gamma,
        0.0,  # to test collisions only, we set the radiative recombination rate to 0
        nlte_raw_plasma_dilution_factor_0_root.alpha_stim,
        nlte_raw_plasma_dilution_factor_0_root.coll_ion_coeff,
        nlte_raw_plasma_dilution_factor_0_root.coll_recomb_coeff,
        nlte_raw_plasma_dilution_factor_0_root.partition_function,
        nlte_raw_plasma_dilution_factor_0_root.levels,
        nlte_raw_plasma_dilution_factor_0_root.level_boltzmann_factor,
        nlte_raw_plasma_dilution_factor_0_root.phi,
        nlte_raw_plasma_dilution_factor_0_root.rate_matrix_index,
        nlte_raw_plasma_dilution_factor_0_root.number_density,
        nlte_raw_plasma_dilution_factor_0_root.nlte_excitation_species,
    )[0]
    ion_number_density_nlte = ion_number_density_nlte.values
    # ion_number_density_nlte[np.abs(ion_number_density_nlte) < 1e-10] = 0.0

    ind = IonNumberDensity(nlte_raw_plasma_dilution_factor_0_root)
    ion_number_density_lte = ind.calculate(
        nlte_raw_plasma_dilution_factor_0_root.thermal_phi_lte,
        nlte_raw_plasma_dilution_factor_0_root.partition_function,
        nlte_raw_plasma_dilution_factor_0_root.number_density,
    )[0]

    ion_number_density_lte = ion_number_density_lte.values
    npt.assert_allclose(
        ion_number_density_lte, ion_number_density_nlte, rtol=1e-2, atol=1e-10
    )


@pytest.mark.xfail
def test_critical_case_dilution_factor_0_lu(
    nlte_raw_plasma_dilution_factor_0_lu,
):
    """Check that the LTE and NLTE solution agree for dilution_factor=0.0."""
    nlte_solver = NLTEPopulationSolverLU(nlte_raw_plasma_dilution_factor_0_lu)
    ion_number_density_nlte = nlte_solver.calculate(
        nlte_raw_plasma_dilution_factor_0_lu.gamma,
        0.0,  # to test collisions only, we set the radiative recombination rate to 0
        nlte_raw_plasma_dilution_factor_0_lu.alpha_stim,
        nlte_raw_plasma_dilution_factor_0_lu.coll_ion_coeff,
        nlte_raw_plasma_dilution_factor_0_lu.coll_recomb_coeff,
        nlte_raw_plasma_dilution_factor_0_lu.partition_function,
        nlte_raw_plasma_dilution_factor_0_lu.levels,
        nlte_raw_plasma_dilution_factor_0_lu.level_boltzmann_factor,
        nlte_raw_plasma_dilution_factor_0_lu.phi,
        nlte_raw_plasma_dilution_factor_0_lu.rate_matrix_index,
        nlte_raw_plasma_dilution_factor_0_lu.number_density,
        nlte_raw_plasma_dilution_factor_0_lu.nlte_excitation_species,
    )[0]
    ion_number_density_nlte = ion_number_density_nlte.values
    # ion_number_density_nlte[np.abs(ion_number_density_nlte) < 1e-10] = 0.0

    ind = IonNumberDensity(nlte_raw_plasma_dilution_factor_0_lu)
    ion_number_density_lte = ind.calculate(
        nlte_raw_plasma_dilution_factor_0_lu.thermal_phi_lte,
        nlte_raw_plasma_dilution_factor_0_lu.partition_function,
        nlte_raw_plasma_dilution_factor_0_lu.number_density,
    )[0]

    ion_number_density_lte = ion_number_density_lte.values
    ion_number_density_lte[
        ion_number_density_lte < 1e-10
    ] = 0.0  # getting rid of small numbers.
    npt.assert_allclose(
        ion_number_density_lte,
        ion_number_density_nlte,
        rtol=1e-2,
        atol=1e-10,
    )
