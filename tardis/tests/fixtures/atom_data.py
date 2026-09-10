from copy import deepcopy

import pytest

from tardis.io.atom_data.base import AtomData

DEFAULT_ATOM_DATA_MD5 = "5d80fa4ae0638469bf1ff281b6ca2a94"


@pytest.fixture(scope="session")
def atomic_data_fname(tardis_regression_path):
    atomic_data_fname = (
        tardis_regression_path
        / "atom_data"
        / "kurucz_cd23_chianti_H_He_latest.h5"
    )

    atom_data_missing_str = (
        f"{atomic_data_fname} atomic datafiles does not seem to exist"
    )

    if not atomic_data_fname.exists():
        pytest.exit(atom_data_missing_str)

    return atomic_data_fname


@pytest.fixture(scope="session")
def atomic_dataset(atomic_data_fname):
    atomic_data = AtomData.from_hdf(atomic_data_fname)

    if atomic_data.md5 != DEFAULT_ATOM_DATA_MD5:
        pytest.skip(
            f'Need default Kurucz atomic dataset (MD5="{DEFAULT_ATOM_DATA_MD5}")'
        )
    else:
        return atomic_data


@pytest.fixture
def kurucz_atomic_data(atomic_dataset):
    atomic_data = deepcopy(atomic_dataset)
    return atomic_data


@pytest.fixture  # (scope="session")
def nlte_atomic_data_fname(tardis_regression_path):
    """
    File name for atomic data used in equilibrium NLTE tests.
    """
    atomic_data_fname = (
        tardis_regression_path
        / "atom_data"
        / "nlte_atom_data"
        / "TestNLTE_He_Ti.h5"
    )

    atom_data_missing_str = (
        f"{atomic_data_fname} atomic datafiles does not seem to exist"
    )

    if not atomic_data_fname.exists():
        pytest.exit(atom_data_missing_str)

    return atomic_data_fname


@pytest.fixture  # (scope="session")
def nlte_atomic_dataset(nlte_atomic_data_fname):
    """
    Atomic dataset used in equilibrium NLTE tests.
    """
    nlte_atomic_data = AtomData.from_hdf(nlte_atomic_data_fname)
    return nlte_atomic_data


@pytest.fixture  # (scope="session")
def nlte_atom_data(nlte_atomic_dataset):
    atomic_data = deepcopy(nlte_atomic_dataset)
    return atomic_data
