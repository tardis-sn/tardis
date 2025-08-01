import logging
import shutil
from pathlib import Path

import yaml
from astropy.config import get_config_dir

from tardis import __path__ as TARDIS_PATH

TARDIS_PATH = Path(TARDIS_PATH[0])
DEFAULT_CONFIG_PATH = (
    TARDIS_PATH / "data" / "default_tardis_internal_config.yml"
)

DEFAULT_DATA_DIR = Path.home() / "Downloads" / "tardis-data"

logger = logging.getLogger(__name__)


def get_internal_configuration():
    config_fpath = Path(get_config_dir()) / "tardis_internal_config.yml"
    if not config_fpath.exists():
        logger.warning(
            f"Configuration File {config_fpath} does not exist - creating new one from default"
        )
        shutil.copy(DEFAULT_CONFIG_PATH, config_fpath)
    with open(config_fpath) as config_fh:
        return yaml.load(config_fh, Loader=yaml.CLoader)


def get_data_dir():
    config = get_internal_configuration()
    data_dir = config.get("data_dir", None)
    if data_dir is None:
        config_fpath = Path(get_config_dir()) / "tardis_internal_config.yml"
        logging.critical(
            f"\n{'*' * 80}\n\nTARDIS will download different kinds of data (e.g. atomic) to its data directory {DEFAULT_DATA_DIR}\n\n"
            f"TARDIS DATA DIRECTORY not specified in {config_fpath}:\n\n"
            f"ASSUMING DEFAULT DATA DIRECTORY {DEFAULT_DATA_DIR}\n "
            f"YOU CAN CHANGE THIS AT ANY TIME IN {config_fpath} \n\n"
            f"{'*' * 80} \n\n"
        )
        if not DEFAULT_DATA_DIR.exists():
            DEFAULT_DATA_DIR.mkdir(parents=True, exist_ok=True)
        config["data_dir"] = DEFAULT_DATA_DIR
        yaml.dump(config, open(config_fpath, "w"), default_flow_style=False)
        data_dir = DEFAULT_DATA_DIR

    if not Path(data_dir).exists():
        raise OSError(f"Data directory specified in {data_dir} does not exist")

    return Path(data_dir)
