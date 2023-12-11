import logging

import yaml

from tardis.io.configuration.config_internal import get_data_dir
from tardis.io.util import download_from_url, get_internal_data_path

logger = logging.getLogger(__name__)


def get_atomic_repo_config():
    """
    Get the repo configuration dictionary for the atomic data

    Returns
    -------
        : dict
    """
    atomic_repo_fname = get_internal_data_path("atomic_data_repo.yml")
    return yaml.load(open(atomic_repo_fname), Loader=yaml.CLoader)


def download_atom_data(atomic_data_name=None, force_download=False):
    """
    Download the atomic data from the repository

    Parameters
    ----------
    atomic_data_name : str
        if None

    force_download : bool
        if True, force download even if file exists

    Returns
    -------
        : None
    """
    atomic_repo = get_atomic_repo_config()

    if atomic_data_name is None:
        atomic_data_name = atomic_repo["default"]

    if atomic_data_name not in atomic_repo:
        raise ValueError(f"Atomic Data name {atomic_data_name} not known")

    dst_fname = get_data_dir() / f"{atomic_data_name}.h5"

    if dst_fname.exists() and not force_download:
        logger.warning(
            f"Atomic Data {atomic_data_name} already exists in {dst_fname}. Will not download - override with force_download=True."
        )
        return
    src_url = atomic_repo[atomic_data_name]["url"]
    mirrors = tuple(atomic_repo[atomic_data_name]["mirrors"])
    checksum = atomic_repo[atomic_data_name]["md5"]

    logger.info(f"Downloading atomic data from {src_url} to {dst_fname}")
    download_from_url(src_url, dst_fname, checksum, mirrors)
