import os
import logging

from pathlib import Path

from tardis.io.configuration.config_internal import get_data_dir
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
    fname : str
        name or path of atom data HDF file

    Returns
    -------
        : str
        resolved fpath
    """
    file_path = Path(fname)

    if file_path.exists():
        return file_path

   # fpath = os.path.join(os.path.join(get_data_dir(), fname))
    data_dir = Path(get_data_dir())
    combined_path = data_dir / file_path
    
    if combined_path.exists():
        logger.info(
            f"\n\tAtom Data {file_path.exist()} not found in local path.\n\tExists in TARDIS Data repo {combined_path}"
        )
        return combined_path

    atom_data_name = file_path.with_suffix("").name
    atom_repo_config = get_atomic_repo_config()
    if atom_data_name in atom_repo_config:
        raise IOError(
            f"Atom Data {fname} not found in path or in TARDIS data repo - it is available as download:\n"
            f"from tardis.io.atom_data.util import download_atom_data\n"
            f"download_atom_data('{atom_data_name}')"
        )

    raise IOError(
        f"Atom Data {fname} is not found in current path or in TARDIS data repo. {atom_data_name} "
        "is also not a standard known TARDIS atom dataset."
    )
