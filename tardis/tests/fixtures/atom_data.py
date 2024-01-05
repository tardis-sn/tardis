from copy import deepcopy
from pathlib import Path
import pytest

from tardis.io.atom_data.base import AtomData
from tardis.io.configuration.config_reader import Configuration
from tardis.model.base import SimulationState

DEFAULT_ATOM_DATA_UUID = "864f1753714343c41f99cb065710cace"


@pytest.fixture(scope="session")
def atomic_data_fname(tardis_ref_path):
    atomic_data_fname = (
        tardis_ref_path / "atom_data" / "kurucz_cd23_chianti_H_He.h5"
    )

    atom_data_missing_str = (
        f"{atomic_data_fname} atomic datafiles " f"does not seem to exist"
    )

    if not atomic_data_fname.exists():
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
    atomic_data_fname = tardis_ref_path / "nlte_atom_data" / "TestNLTE_He_Ti.h5"

    atom_data_missing_str = (
        f"{atomic_data_fname} atomic datafiles " f"does not seem to exist"
    )

    if not atomic_data_fname.exists():
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


@pytest.fixture  # (scope="session")
def tardis_model_config_nlte_root(example_configuration_dir):
    config = Configuration.from_yaml(
        example_configuration_dir / "tardis_configv1_nlte.yml"
    )
    config.plasma.nlte_solver = "root"
    return config


@pytest.fixture  # (scope="session")
def tardis_model_config_nlte_lu(example_configuration_dir):
    config = Configuration.from_yaml(
        example_configuration_dir / "tardis_configv1_nlte.yml"
    )
    config.plasma.nlte_solver = "lu"
    return config


@pytest.fixture  # (scope="session")
def nlte_raw_model_root(tardis_model_config_nlte_root, nlte_atom_data):
    return SimulationState.from_config(
        tardis_model_config_nlte_root, nlte_atom_data
    )


@pytest.fixture  # (scope="session")
def nlte_raw_model_lu(tardis_model_config_nlte_lu, nlte_atom_data):
    return SimulationState.from_config(
        tardis_model_config_nlte_lu, nlte_atom_data
    )
