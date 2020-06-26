import os
import logging

from tardis.io.config_internal import get_data_dir
from tardis.io.atom_data.atom_web_download import (
    get_atomic_repo_config,
    download_atom_data,
)

logger = logging.getLogger(__name__)


def resolve_atom_data_fname(fname):
    """
    Check where if atom data HDF file is available on disk, can be downloaded or does not exist

    Parameters
    ----------
    fname: str
        name or path of atom data HDF file

    Returns
    -------
        : str
        resolved fpath
    """

    if os.path.exists(fname):
        return fname

    fpath = os.path.join(os.path.join(get_data_dir(), fname))
    if os.path.exists(fpath):
        logger.info(
            "Atom Data {0} not found in local path. Exists in TARDIS Data repo {1}".format(
                fname, fpath
            )
        )
        return fpath

    atom_data_name = fname.replace(".h5", "")
    atom_repo_config = get_atomic_repo_config()
    if atom_data_name in atom_repo_config:
        raise IOError(
            "Atom Data {0} not found in path or in TARDIS data repo - it is available as download:\n"
            "from tardis.io.atom_data import download_atom_data\n"
            "download_atom_data('{1}')".format(fname, atom_data_name)
        )

    raise IOError(
        f"Atom Data {fname} is not found in current path or in TARDIS data repo. {atom_data_name} "
        "is also not a standard known TARDIS atom dataset."
    )
