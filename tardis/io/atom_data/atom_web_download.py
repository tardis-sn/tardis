import os
import logging

from tardis.io.util import get_internal_data_path, download_from_url
from tardis.io.config_internal import get_data_dir
import yaml

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


def download_atom_data(atomic_data_name=None):
    """
    Download the atomic data from the repository

    Parameters
    ----------
    atomic_data_name: str
        if None

    Returns
    -------
        : None

    """
    atomic_repo = get_atomic_repo_config()

    if atomic_data_name is None:
        atomic_data_name = atomic_repo["default"]

    if atomic_data_name not in atomic_repo:
        raise ValueError(
            "Atomic Data name {0} not known".format(atomic_data_name)
        )
    dst_dir = os.path.join(get_data_dir(), "{0}.h5".format(atomic_data_name))
    src_url = atomic_repo[atomic_data_name]["url"]
    logger.info(
        "Downloading atomic data from {0} to {1}".format(src_url, dst_dir)
    )
    download_from_url(src_url, dst_dir)
