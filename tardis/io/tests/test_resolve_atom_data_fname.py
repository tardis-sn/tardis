import pytest
import os

from tardis.io.atom_data.util import resolve_atom_data_fname


@pytest.fixture
def atom_data_fname_absolute(tmp_path):
    """
    Creates a temp absolute path to mock fname.
    """
    fname = tmp_path / "mock_atom_data.h5"
    fname.touch()
    return fname


def test_absolute_path_atom_data(atom_data_fname_absolute):
    """
    Checks whether the absolute path provided in the config exists or not.

    Parameters
    ----------
    atom_data_fname_absolute : tmp_path
        Absolute path to the atomic data
    """
    fpath = str(atom_data_fname_absolute)
    assert resolve_atom_data_fname(fpath) == fpath

    non_existent_fpath = str(
        atom_data_fname_absolute.with_name("non_existent_file.txt")
    )
    with pytest.raises(IOError):
        resolve_atom_data_fname(non_existent_fpath)


@pytest.fixture
def data_path(tmp_path):
    """
    Creates a temporary directory that contains the atomic data.
    """
    # Create a temporary directory
    data_dir = tmp_path / "data_dir"
    data_dir.mkdir()
    fname = data_dir / "mock_atom_data.h5"
    fname.touch()
    return data_dir


def test_relative_atom_data_fname(data_path, monkeypatch):
    """
    Checks whether the file name in config is searched in the atomic-data directory.

    Parameters
    ----------
    data_path : tmp_path
        Atomic Data directory that is used to stub the return value of 'get_data_dir()'
    monkeypatch : Monkeypatch pytest fixture
        Replaces get_data_dir() with a mock function
    """

    monkeypatch.setattr(
        "tardis.io.atom_data.util.get_data_dir", lambda: str(data_path)
    )
    resolved_fname = resolve_atom_data_fname("mock_atom_data.h5")

    expected_fname = os.path.join(str(data_path), "mock_atom_data.h5")
    assert resolved_fname == expected_fname

    with pytest.raises(IOError):
        resolve_atom_data_fname("non_existent_data.h5")
