from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
import pytest

from tardis.iip_plasma.continuum.base_continuum_data import ContinuumData
from tardis.iip_plasma.standard_plasmas import LegacyPlasmaArray
from tardis.io.atom_data import AtomData
from tardis.io.configuration.config_reader import Configuration


@pytest.fixture
def iip_atom_data(tardis_regression_path):
    # identical atomic data to that used by C Vogl
    atom_data = AtomData.from_hdf(
        tardis_regression_path
        / "atom_data"
        / "christians_atomdata_converted_04Dec25.h5"
    )  # currently not available for public use

    atom_data.prepare_atom_data([1], "macroatom", [(1, 0)], [(1, 0)])

    atom_data.continuum_data = ContinuumData(
        atom_data, selected_continuum_species=[(1, 0)]
    )

    atom_data.yg_data.columns = list(atom_data.collision_data_temperatures)

    atom_data.nlte_data._init_indices()

    atom_data.has_collision_data = False

    return atom_data


@pytest.fixture
def ctardis_compare_config():
    config = Configuration.from_yaml(
        "tardis/workflows/tests/data/ctardis_compare.yml"
    )
    return config


@pytest.fixture
def elemental_number_density(tardis_regression_path):
    elemental_number_density = pd.read_hdf(
        tardis_regression_path
        / "tardis"
        / "workflows"
        / "tests"
        / "ctardis_elemental_number_density.h5",
        key="data",
    )
    elemental_number_density.columns = elemental_number_density.columns.astype(
        int
    )
    return elemental_number_density


@pytest.fixture
def iip_plasma(iip_atom_data, elemental_number_density, ctardis_compare_config):
    plasma = LegacyPlasmaArray(
        elemental_number_density,
        iip_atom_data,
        ctardis_compare_config.supernova.time_explosion.to("s").value,
        nlte_config=ctardis_compare_config.plasma.nlte,
        delta_treatment=None,
        ionization_mode="nlte",
        excitation_mode="dilute-lte",
        line_interaction_type=ctardis_compare_config.plasma.line_interaction_type,
        link_t_rad_t_electron=1.0 * np.ones(24),
        # link_t_rad_t_electron=self.ws**0.25,
        helium_treatment="none",
        heating_rate_data_file=None,
        v_inner=None,
        v_outer=None,
        continuum_treatment=True,
    )

    return plasma


def test_iip_plasma_initialization(
    iip_plasma, ctardis_compare_config, tardis_regression_path
):
    j_blues_ctardis = pd.read_hdf(
        tardis_regression_path
        / "tardis"
        / "workflows"
        / "tests"
        / "ctardis_j_blues_ctardis_init_nlte.h5",
        key="data",
    )

    radiation_temp = 9984.96131287 * np.ones(24)
    dilution_factor = np.array(
        [
            0.18635244,
            0.15938095,
            0.11736085,
            0.34665656,
            0.32265696,
            0.30224056,
            0.28436446,
            0.26841929,
            0.2540108,
            0.24086562,
            0.22878441,
            0.21761613,
            0.20724285,
            0.1975702,
            0.18852112,
            0.18003167,
            0.17204798,
            0.16452412,
            0.15742053,
            0.15070279,
            0.14434073,
            0.13830767,
            0.13257993,
            0.12856901,
        ]
    )

    iip_plasma.update_radiationfield(
        radiation_temp,
        dilution_factor,
        j_blues_ctardis,
        ctardis_compare_config.plasma.nlte,
        initialize_nlte=True,
        n_e_convergence_threshold=0.05,
        **{},
    )

    tau_sobolevs_ctardis = pd.read_hdf(
        tardis_regression_path
        / "tardis"
        / "workflows"
        / "tests"
        / "ctardis_tau_sobolevs_init_nlte.h5",
        key="data",
    )

    beta_sobolevs_ctardis = pd.read_hdf(
        tardis_regression_path
        / "tardis"
        / "workflows"
        / "tests"
        / "ctardis_beta_sobolevs_init_nlte.h5",
        key="data",
    )
    ion_number_density_ctardis = pd.read_hdf(
        tardis_regression_path
        / "tardis"
        / "workflows"
        / "tests"
        / "ctardis_ion_density_init_nlte.h5",
        key="data",
    )
    level_number_density_ctardis = pd.read_hdf(
        tardis_regression_path
        / "tardis"
        / "workflows"
        / "tests"
        / "ctardis_level_number_density_init_nlte.h5",
        key="data",
    )
    transition_probabilities_ctardis = pd.read_hdf(
        tardis_regression_path
        / "tardis"
        / "workflows"
        / "tests"
        / "ctardis_transition_probabilities_init_nlte.h5",
        key="data",
    )

    pd.testing.assert_frame_equal(
        iip_plasma.transition_probabilities,
        transition_probabilities_ctardis,
        rtol=1e-7,
        atol=0,
        check_dtype=False,
    )

    pd.testing.assert_frame_equal(
        iip_plasma.ion_number_density,
        ion_number_density_ctardis,
        rtol=1e-6,
        atol=0,
        check_dtype=False,
    )

    # Sobolev values are stored differently between codes, so comparing raw data instead
    np.testing.assert_allclose(
        iip_plasma.tau_sobolevs.values,
        tau_sobolevs_ctardis.values,
        rtol=1e-6,
        atol=0,
    )

    np.testing.assert_allclose(
        iip_plasma.beta_sobolev.values,
        beta_sobolevs_ctardis.values,
        rtol=1e-6,
        atol=0,
    )

    pd.testing.assert_frame_equal(
        iip_plasma.level_number_density,
        level_number_density_ctardis,
        rtol=1e-7,
        atol=0,
        check_dtype=False,
    )
