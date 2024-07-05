from pathlib import Path
import numpy as np
import pandas as pd
import numpy.testing as npt

from astropy import units as u
from tardis.io.configuration.config_reader import Configuration
from tardis.io.atom_data.base import AtomData
from tardis.model import SimulationState
import pytest


@pytest.fixture(
    scope="module",
    params=[
        "model_full",
        "branch85",
        "uniform",
        "powerlaw",
        "exponential",
        "radiative",
    ],
)
def model_config_fnames(request, example_csvy_file_dir):
    """Function to retrieve filenames of target data for tests"""
    csvy_config_file = example_csvy_file_dir / f"{request.param}_csvy.yml"
    old_config_file = example_csvy_file_dir / f"{request.param}_old_config.yml"
    return csvy_config_file, old_config_file


def test_compare_models(model_config_fnames, atomic_dataset):
    """Compare identical models produced by .from_config and
    .from_csvy to check that velocities, densities and abundances
    (pre and post decay) are the same"""
    csvy_config_file, old_config_file = model_config_fnames
    tardis_config = Configuration.from_yaml(csvy_config_file)
    tardis_config_old = Configuration.from_yaml(old_config_file)
    csvy_simulation_state = SimulationState.from_csvy(
        tardis_config, atom_data=atomic_dataset
    )
    config_simulation_state = SimulationState.from_config(
        tardis_config_old, atom_data=atomic_dataset
    )
    csvy_model_props = csvy_simulation_state.get_properties().keys()
    config_model_props = config_simulation_state.get_properties().keys()
    npt.assert_array_equal(csvy_model_props, config_model_props)
    for prop in config_model_props:
        csvy_model_val = csvy_simulation_state.get_properties()[prop]
        config_model_val = config_simulation_state.get_properties()[prop]
        if prop == "homologous_density":
            npt.assert_array_almost_equal(
                csvy_model_val.density_0.value, config_model_val.density_0.value
            )
            npt.assert_array_almost_equal(
                csvy_model_val.time_0.value, config_model_val.time_0.value
            )
        else:
            if hasattr(config_model_val, "value"):
                config_model_val = config_model_val.value
                csvy_model_val = csvy_model_val.value
            npt.assert_array_almost_equal(csvy_model_val, config_model_val)

    assert (
        csvy_simulation_state.abundance.shape
        == config_simulation_state.abundance.shape
    )
    assert (
        csvy_simulation_state.composition.nuclide_mass_fraction.shape
        == config_simulation_state.composition.nuclide_mass_fraction.shape
    )
    assert (
        csvy_simulation_state.abundance.shape
        == config_simulation_state.abundance.shape
    )
    npt.assert_array_almost_equal(
        csvy_simulation_state.abundance.to_numpy(),
        config_simulation_state.abundance.to_numpy(),
    )
    npt.assert_array_almost_equal(
        csvy_simulation_state.composition.nuclide_mass_fraction.to_numpy(),
        config_simulation_state.composition.nuclide_mass_fraction.to_numpy(),
    )
    npt.assert_array_almost_equal(
        csvy_simulation_state.abundance.to_numpy(),
        config_simulation_state.abundance.to_numpy(),
    )


def test_dimensionality_after_update_v_inner_boundary(
    example_csvy_file_dir, atomic_dataset
):
    """Test that the dimensionality of SimulationState parameters after updating v_inner_boundary
    in the context of csvy models specifically"""
    csvy_config_file = example_csvy_file_dir / "radiative_csvy.yml"
    config = Configuration.from_yaml(csvy_config_file)
    csvy_model = SimulationState.from_csvy(config, atom_data=atomic_dataset)

    new_config = config
    new_config.model.v_inner_boundary = csvy_model.velocity[1]
    new_csvy_model = SimulationState.from_csvy(
        new_config, atom_data=atomic_dataset
    )

    assert new_csvy_model.no_of_raw_shells == csvy_model.no_of_raw_shells
    assert new_csvy_model.no_of_shells == csvy_model.no_of_shells - 1
    assert new_csvy_model.velocity.shape[0] == csvy_model.velocity.shape[0] - 1
    assert new_csvy_model.density.shape[0] == csvy_model.density.shape[0] - 1
    assert new_csvy_model.volume.shape[0] == csvy_model.volume.shape[0] - 1
    assert (
        new_csvy_model.t_radiative.shape[0]
        == csvy_model.t_radiative.shape[0] - 1
    )


@pytest.fixture(scope="module")
def csvy_model_test_abundances(example_csvy_file_dir, atomic_dataset):
    """Returns SimulationState to use to test abundances dataframes"""
    csvypath = example_csvy_file_dir / "csvy_model_to_test_abundances.yml"
    config = Configuration.from_yaml(csvypath)
    csvy_model_test_abundances = SimulationState.from_csvy(
        config, atom_data=atomic_dataset
    )
    return csvy_model_test_abundances


@pytest.fixture(scope="function")
def reference_input_dataframes():
    """Reference elemental and isotope abundances dataframes for
    `test_read_csvy`.
    Constructed from `csvy_model_to_test_abundances.yml`
    Rows in dataframes are abundances for a fixed element or isotope;
    columns represent different shells"""

    abundance_index = pd.Index([1, 2], name="atomic_number")
    reference_input_abundance = pd.DataFrame(
        [[0.0, 0.33, 0.3, 0.5, 0.4, 0.2], [0.98, 0.64, 0.6, 0.4, 0.55, 0.79]],
        index=abundance_index,
    )

    arrays = [[28], [56]]
    isotope_index = pd.MultiIndex.from_arrays(
        arrays, names=["atomic_number", "mass_number"]
    )
    reference_input_isotopes = pd.DataFrame(
        [[0.02, 0.03, 0.1, 0.1, 0.05, 0.01]],
        columns=np.arange(6),
        index=isotope_index,
    )
    return reference_input_abundance, reference_input_isotopes


def test_read_csvy_abundances(
    csvy_model_test_abundances, reference_input_dataframes
):
    """Test if model reads abundances and isotope abundances
    and constructs dataframes correctly before applying decay"""
    (
        reference_input_abundance,
        reference_input_isotopes,
    ) = reference_input_dataframes

    composition = csvy_model_test_abundances.composition
    nuclide_mass_fraction = composition.nuclide_mass_fraction
    model_abundances = nuclide_mass_fraction[
        nuclide_mass_fraction.index.get_level_values(1) == -1
    ]

    reference_input_shape = reference_input_abundance.shape
    assert model_abundances.shape == reference_input_shape
    npt.assert_array_almost_equal(
        reference_input_abundance.to_numpy(),
        model_abundances.to_numpy(),
    )

    model_isotopes = composition.isotopic_mass_fraction
    # reference_input_isotopes_shape = reference_input_isotopes.shape
    # We can't assert the shape anymore because the isotope abundances used
    # to be decayed after being loaded into SimulationState.
    #  Now the abundances are decayed before
    # assert model_isotopes.shape == reference_input_isotopes_shape
    # Same applies to the comparison - we are summing up the mass_fractions to compare pre/post decay
    npt.assert_array_almost_equal(
        reference_input_isotopes.to_numpy()[0],
        model_isotopes.sum(axis=0).to_numpy(),
        decimal=1,
    )


@pytest.fixture(scope="function")
def reference_decayed_abundance():
    """Returns reference abundances dataframe for `test_csvy_model_decay`
    Constructed from `csvy_model_to_test_abundances.yml`
    For the decay calculations the following procedure is used:
    Ni_halflife = 6.075 * u.d
    Co_halflife = 77.233 * u.d

    lambda_Ni = np.log(2) / Ni_halflife
    lambda_Co = np.log(2) / Co_halflife

    t = 4 * u.d means 4 days have passed since the explosion

    def N1(N0, lambda1, t=4.0 * u.d):
        return N0 * np.exp(-lambda1 * t)


    def N2(N1_0, lambda_1, lambda_2, t=4.0 * u.d):
        return (lambda_1 * N1_0
                * (np.exp(-lambda_1 * t) / (lambda_2 - lambda_1)
                + np.exp(-lambda_2 * t) / (lambda_1 - lambda_2)))

     if the original Ni56 abundance for a given shell is 0.05, after 4 days:

     cobalt_abundance_after_4_days = N2(0.05, lambda_Ni, lambda_Co)
     nickel_abundance_after_4_days = N1(0.05, lambda_Ni)
     iron_abundance_after_4_days = 0.05 - cobalt_abundance_after_4_days
                                    - nickel_abundance_after_4_days
     In the reference_decayed_dataframe every row represents a specific element
     and every column represents a shell"""
    decay_index = pd.Index([1, 2, 26, 27, 28], name="atomic_number")
    reference_decayed_abundance = pd.DataFrame(
        [
            [0.0, 0.33, 0.3, 0.5, 0.4, 0.2],
            [0.98, 0.64, 0.6, 0.4, 0.55, 0.79],
            [
                0.00013977843354947162,
                0.00020966765032420787,
                0.0006988921677473642,
                0.0006988921677473642,
                0.0003494460838736821,
                6.988921677473581e-05,
            ],
            [
                0.007188928223953217,
                0.010783392335929825,
                0.035944641119766085,
                0.035944641119766085,
                0.017972320559883043,
                0.0035944641119766084,
            ],
            [
                0.012671293342497312,
                0.019006940013745966,
                0.06335646671248656,
                0.06335646671248656,
                0.03167823335624328,
                0.006335646671248656,
            ],
        ],
        index=decay_index,
    )
    return reference_decayed_abundance


def test_csvy_model_decay(
    csvy_model_test_abundances, reference_decayed_abundance
):
    """Compare model abundance decay against decay calculations
    done by hand."""
    model_decayed_abundance_shape = csvy_model_test_abundances.abundance.shape
    reference_decayed_abundance_shape = reference_decayed_abundance.shape
    assert model_decayed_abundance_shape == reference_decayed_abundance_shape
    npt.assert_array_almost_equal(
        reference_decayed_abundance.to_numpy(),
        csvy_model_test_abundances.abundance.to_numpy(),
    )
