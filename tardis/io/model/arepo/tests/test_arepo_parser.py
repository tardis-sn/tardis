import json
from pathlib import Path

import numpy as np
import pytest
from astropy import units as u

from tardis.io.model.arepo import ArepoData
from tardis.io.model.arepo.utils import (
    create_cone_profile,
    create_full_profile,
    export_profile_to_csvy,
    rebin_profile,
)


@pytest.fixture
def arepo_snapshot_fname(tardis_regression_path):
    return Path(tardis_regression_path) / "arepo_data" / "arepo_snapshot.json"


@pytest.fixture
def get_cone_csvy_model(arepo_snapshot_fname, tmp_path):
    with open(arepo_snapshot_fname) as json_file:
        snapshot_data = json.loads(json.load(json_file))

    position, velocities, densities, mass, isotope_fractions, time = (
        snapshot_data["pos"],
        snapshot_data["vel"],
        snapshot_data["rho"],
        snapshot_data["mass"],
        snapshot_data["xnuc"],
        snapshot_data["time"],
    )
    position = np.array(position)
    velocities = np.array(velocities)
    densities = np.array(densities)
    mass = np.array(mass)

    # The nuclear data should be in a dict where each element has its own entry (with the key being the element name)
    isotope_dict = {
        "ni56": np.array(isotope_fractions[0]),
        "si28": np.array(isotope_fractions[1]),
    }

    arepo_data = ArepoData(
        time=u.Quantity(time, "s"),
        position=position,
        velocities=velocities,
        densities=densities,
        mass=mass,
        isotope_dict=isotope_dict,
    )

    (
        pos_prof_p,
        pos_prof_n,
        vel_prof_p,
        vel_prof_n,
        rho_prof_p,
        rho_prof_n,
        mass_prof_p,
        mass_prof_n,
        xnuc_prof_p,
        xnuc_prof_n,
    ) = create_cone_profile(
        arepo_data, opening_angle=40, inner_radius=1e11, outer_radius=2e11
    )

    (
        pos_prof_p,
        pos_prof_n,
        vel_prof_p,
        vel_prof_n,
        rho_prof_p,
        rho_prof_n,
        mass_prof_p,
        mass_prof_n,
        xnuc_prof_p,
        xnuc_prof_n,
    ) = rebin_profile(
        pos_prof_p,
        pos_prof_n,
        vel_prof_p,
        vel_prof_n,
        rho_prof_p,
        rho_prof_n,
        mass_prof_p,
        mass_prof_n,
        xnuc_prof_p,
        xnuc_prof_n,
        nshells=20,
    )

    testfile = export_profile_to_csvy(
        pos_prof_p,
        vel_prof_p,
        rho_prof_p,
        mass_prof_p,
        xnuc_prof_p,
        arepo_data.time,
        tmp_path / "arepo_parser_test.csvy",
        nshells=20,
    )

    with open(testfile) as file:
        csvy_lines = file.readlines()

    return csvy_lines


@pytest.fixture
def get_full_csvy_model(arepo_snapshot_fname, tmp_path):
    with open(arepo_snapshot_fname) as json_file:
        snapshot_data = json.loads(json.load(json_file))

    position, velocities, densities, mass, isotope_fractions, time = (
        snapshot_data["pos"],
        snapshot_data["vel"],
        snapshot_data["rho"],
        snapshot_data["mass"],
        snapshot_data["xnuc"],
        snapshot_data["time"],
    )
    position = np.array(position)
    velocities = np.array(velocities)
    densities = np.array(densities)
    mass = np.array(mass)

    # The nuclear data should be in a dict where each element has its own entry (with the key being the element name)
    isotope_dict = {
        "ni56": np.array(isotope_fractions[0]),
        "si28": np.array(isotope_fractions[1]),
    }

    arepo_data = ArepoData(
        time=u.Quantity(time, "s"),
        position=position,
        velocities=velocities,
        densities=densities,
        mass=mass,
        isotope_dict=isotope_dict,
    )

    (
        pos_prof_p,
        pos_prof_n,
        vel_prof_p,
        vel_prof_n,
        rho_prof_p,
        rho_prof_n,
        mass_prof_p,
        mass_prof_n,
        xnuc_prof_p,
        xnuc_prof_n,
    ) = create_full_profile(arepo_data, inner_radius=1e11, outer_radius=2e11)

    (
        pos_prof_p,
        pos_prof_n,
        vel_prof_p,
        vel_prof_n,
        rho_prof_p,
        rho_prof_n,
        mass_prof_p,
        mass_prof_n,
        xnuc_prof_p,
        xnuc_prof_n,
    ) = rebin_profile(
        pos_prof_p,
        pos_prof_n,
        vel_prof_p,
        vel_prof_n,
        rho_prof_p,
        rho_prof_n,
        mass_prof_p,
        mass_prof_n,
        xnuc_prof_p,
        xnuc_prof_n,
        nshells=20,
    )

    testfile = export_profile_to_csvy(
        pos_prof_p,
        vel_prof_p,
        rho_prof_p,
        mass_prof_p,
        xnuc_prof_p,
        arepo_data.time,
        tmp_path / "arepo_parser_test.csvy",
        nshells=20,
    )

    with open(testfile) as file:
        csvy_lines = file.readlines()

    return csvy_lines


def test_cone_profile(get_cone_csvy_model, regression_data):
    """
    Test cone profile generation against regression data.

    Parameters
    ----------
    get_cone_csvy_model : list
        Generated CSVY lines from cone profile.
    regression_data : pytest fixture
        Regression data fixture for comparison.
    """
    csvy_array = np.array(
        get_cone_csvy_model
    )  # this is essentially a file (list of lines)
    expected = regression_data.sync_ndarray(csvy_array)
    np.testing.assert_array_equal(csvy_array, expected)


def test_full_profile(get_full_csvy_model, regression_data):
    """
    Test full profile generation against regression data.

    Parameters
    ----------
    get_full_csvy_model : list
        Generated CSVY lines from full profile.
    regression_data : pytest fixture
        Regression data fixture for comparison.
    """
    csvy_array = np.array(
        get_full_csvy_model
    )  # this is essentially a file (list of lines)
    expected = regression_data.sync_ndarray(csvy_array)
    np.testing.assert_array_equal(csvy_array, expected)
