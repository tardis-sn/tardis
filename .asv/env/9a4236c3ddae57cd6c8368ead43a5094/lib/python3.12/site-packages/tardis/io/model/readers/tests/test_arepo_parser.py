from pathlib import Path
from tardis.io.model.readers import arepo
import pytest
import numpy as np
import os
import json


@pytest.fixture()
def arepo_snapshot_fname(tardis_ref_path):
    return Path(tardis_ref_path) / "arepo_data" / "arepo_snapshot.json"


@pytest.fixture
def get_cone_csvy_model(arepo_snapshot_fname, example_model_file_dir):
    with open(arepo_snapshot_fname, "r") as json_file:
        data = json.loads(json.load(json_file))

    pos, vel, rho, mass, nucs, time = (
        data["pos"],
        data["vel"],
        data["rho"],
        data["mass"],
        data["xnuc"],
        data["time"],
    )
    pos = np.array(pos)
    vel = np.array(vel)
    rho = np.array(rho)
    mass = np.array(mass)

    # The nuclear data should be in a dict where each element has its own entry (with the key being the element name)
    xnuc = {
        "ni56": np.array(nucs[0]),
        "si28": np.array(nucs[1]),
    }

    profile = arepo.ConeProfile(pos, vel, rho, mass, xnuc, time)

    profile.create_profile(
        opening_angle=40, inner_radius=1e11, outer_radius=2e11
    )

    testfile = profile.export(
        20, example_model_file_dir / "arepo_parser_test.csvy"
    )

    with open(testfile, "r") as file:
        data = file.readlines()

    print(testfile)
    os.remove(testfile)

    return data


@pytest.fixture
def get_full_csvy_model(arepo_snapshot_fname, example_model_file_dir):
    with open(arepo_snapshot_fname, "r") as json_file:
        data = json.loads(json.load(json_file))

    pos, vel, rho, mass, nucs, time = (
        data["pos"],
        data["vel"],
        data["rho"],
        data["mass"],
        data["xnuc"],
        data["time"],
    )
    pos = np.array(pos)
    vel = np.array(vel)
    rho = np.array(rho)
    mass = np.array(mass)

    # The nuclear data should be in a dict where each element has its own entry (with the key being the element name)
    xnuc = {
        "ni56": np.array(nucs[0]),
        "si28": np.array(nucs[1]),
    }

    profile = arepo.FullProfile(pos, vel, rho, mass, xnuc, time)

    profile.create_profile(inner_radius=1e11, outer_radius=2e11)

    testfile = profile.export(
        20, example_model_file_dir / "arepo_parser_test.csvy"
    )

    with open(testfile, "r") as file:
        data = file.readlines()

    os.remove(testfile)

    return data


@pytest.fixture
def get_cone_reference_data(example_model_file_dir):
    with open(
        example_model_file_dir / "arepo_cone_reference_model.csvy", "r"
    ) as file:
        data = file.readlines()

    return data


@pytest.fixture
def get_full_reference_data(example_model_file_dir):
    with open(
        example_model_file_dir / "arepo_full_reference_model.csvy", "r"
    ) as file:
        data = file.readlines()

    return data


@pytest.mark.ignore_generate
def test_cone_profile(get_cone_csvy_model, get_cone_reference_data):
    assert get_cone_csvy_model == get_cone_reference_data


@pytest.mark.ignore_generate
def test_full_profile(get_full_csvy_model, get_full_reference_data):
    assert get_full_csvy_model == get_full_reference_data
