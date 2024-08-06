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
    fname : Path
        name or path of atom data HDF file

    Returns
    -------
        : Path
        resolved fpath
    """

    fname = Path(fname)
    if os.path.exists(fname):
        return fname

    fname = Path(fname.stem).with_suffix(".h5")
    fpath = Path(os.path.join(get_data_dir(), fname))
    if os.path.exists(fpath):
        logger.info(
            f"\n\tAtom Data {fname} not found in local path.\n\tExists in TARDIS Data repo {fpath}"
        )
        return fpath

    atom_data_name = fname.stem
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


def set_atom_data_attributes(atom_data, store, attribute):
    """Sets arbitrary atom data attributes, throws error and sets to None
    if they are not available.

    Parameters
    ----------
    atom_data : AtomData
        The atom data to modify
    store : pd.HDFStore
        Data source
    property : str
        Property to modify
    """
    try:
        setattr(atom_data, attribute, store.root._v_attrs[attribute])
        new_attribute = getattr(atom_data, attribute)
        if hasattr(new_attribute, "decode"):
            setattr(
                atom_data,
                attribute,
                store.root._v_attrs[attribute].decode("ascii"),
            )
    except KeyError:
        logger.debug(
            f"{attribute} not available for Atom Data. Setting value to None"
        )
        setattr(atom_data, attribute, None)
