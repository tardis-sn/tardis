import os
from copy import deepcopy

import pytest

from tardis.io.atom_data.base import AtomData
from tardis.io.config_reader import Configuration
from tardis.model.base import Radial1DModel

DEFAULT_ATOM_DATA_UUID = "864f1753714343c41f99cb065710cace"


@pytest.fixture(scope="session")
def atomic_data_fname(tardis_ref_path):
    atomic_data_fname = os.path.join(
        tardis_ref_path, "atom_data", "kurucz_cd23_chianti_H_He.h5"
    )

    atom_data_missing_str = (
        f"{atomic_data_fname} atomic datafiles " f"does not seem to exist"
    )

    if not os.path.exists(atomic_data_fname):
        pytest.exit(atom_data_missing_str)

    return atomic_data_fname


@pytest.fixture(scope="session")
def atomic_dataset(atomic_data_fname):
    atomic_data = AtomData.from_hdf(atomic_data_fname)

    if atomic_data.md5 != DEFAULT_ATOM_DATA_UUID:
        pytest.skip(
            f'Need default Kurucz atomic dataset (md5="{DEFAULT_ATOM_DATA_UUID}")'
        )
    else:
        return atomic_data


@pytest.fixture
def kurucz_atomic_data(atomic_dataset):
    atomic_data = deepcopy(atomic_dataset)
    return atomic_data


@pytest.fixture  # (scope="session")
def nlte_atomic_data_fname(tardis_ref_path):
    """
    File name for the atomic data file used in NTLE ionization solver tests.
    """
    atomic_data_fname = os.path.join(
        tardis_ref_path, "nlte_atom_data", "TestNLTE_He_Ti.h5"
    )

    atom_data_missing_str = (
        f"{atomic_data_fname} atomic datafiles " f"does not seem to exist"
    )

    if not os.path.exists(atomic_data_fname):
        pytest.exit(atom_data_missing_str)

    return atomic_data_fname


@pytest.fixture  # (scope="session")
def nlte_atomic_dataset(nlte_atomic_data_fname):
    """
    Atomic dataset used for NLTE ionization solver tests.
    """
    nlte_atomic_data = AtomData.from_hdf(nlte_atomic_data_fname)
    return nlte_atomic_data


@pytest.fixture  # (scope="session")
def nlte_atom_data(nlte_atomic_dataset):

    atomic_data = deepcopy(nlte_atomic_dataset)
    return atomic_data


data_path = os.path.join("tardis", "io", "tests", "data")


@pytest.fixture  # (scope="session")
def tardis_model_config_nlte():
    filename = "tardis_configv1_nlte.yml"
    return Configuration.from_yaml(os.path.join(data_path, filename))


@pytest.fixture  # (scope="session")
def nlte_raw_model(tardis_model_config_nlte):
    return Radial1DModel.from_config(tardis_model_config_nlte)
