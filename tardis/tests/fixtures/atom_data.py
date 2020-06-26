import os
from copy import deepcopy

import pytest

from tardis.io.atom_data.base import AtomData

DEFAULT_ATOM_DATA_UUID = "864f1753714343c41f99cb065710cace"


@pytest.fixture(scope="session")
def atomic_data_fname(tardis_ref_path):
    atomic_data_fname = os.path.join(
        tardis_ref_path, "atom_data", "kurucz_cd23_chianti_H_He.h5"
    )

    atom_data_missing_str = (
        "{0} atomic datafiles "
        "does not seem to exist".format(atomic_data_fname)
    )

    if not os.path.exists(atomic_data_fname):
        pytest.exit(atom_data_missing_str)

    return atomic_data_fname


@pytest.fixture(scope="session")
def atomic_dataset(atomic_data_fname):
    atomic_data = AtomData.from_hdf(atomic_data_fname)

    if atomic_data.md5 != DEFAULT_ATOM_DATA_UUID:
        pytest.skip(
            'Need default Kurucz atomic dataset (md5="{}"'.format(
                DEFAULT_ATOM_DATA_UUID
            )
        )
    else:
        return atomic_data


@pytest.fixture
def kurucz_atomic_data(atomic_dataset):
    atomic_data = deepcopy(atomic_dataset)
    return atomic_data
