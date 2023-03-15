import pytest
import os

from tardis.io.atom_data.util import resolve_atom_data_fname


@pytest.fixture
def atom_data_fname_absolute(tmp_path):
    # create a temp file to mock fname
    fname = tmp_path / "mock_atom_data.h5"
    fname.touch()
    return fname


def test_absolute_path_atom_data(atom_data_fname_absolute):
    fpath = str(atom_data_fname_absolute)
    assert resolve_atom_data_fname(fpath) == fpath

    non_existent_fpath = str(
        atom_data_fname_absolute.with_name("non_existent_file.txt")
    )
    with pytest.raises(IOError):
        resolve_atom_data_fname(non_existent_fpath)
