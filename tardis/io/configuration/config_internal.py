import logging
import shutil
from importlib import resources
from pathlib import Path

import yaml
from astropy.config import get_config_dir

TARDIS_PATH = Path(str(resources.files("tardis")))

DEFAULT_CONFIG_PATH = (
    TARDIS_PATH / "data" / "default_tardis_internal_config.yml"
)

DEFAULT_DATA_DIR = Path.home() / "Downloads" / "tardis-data"

logger = logging.getLogger(__name__)


def get_internal_configuration():
    config_fpath = Path(get_config_dir()) / "tardis_internal_config.yml"
    if not config_fpath.exists():
        logger.warning(
            "Configuration File %s does not exist - creating new one from default", 
            config_fpath
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
            "\n%s\n\nTARDIS will download different kinds of data (e.g. atomic) to its data directory %s\n\n"
            "TARDIS DATA DIRECTORY not specified in %s:\n\n"
            "ASSUMING DEFAULT DATA DIRECTORY %s\n "
            "YOU CAN CHANGE THIS AT ANY TIME IN %s \n\n"
            "%s \n\n",
            "*" * 80,
            DEFAULT_DATA_DIR,
            config_fpath,
            DEFAULT_DATA_DIR,
            config_fpath,
            "*" * 80
        )
        if not DEFAULT_DATA_DIR.exists():
            DEFAULT_DATA_DIR.mkdir(parents=True, exist_ok=True)
        config["data_dir"] = DEFAULT_DATA_DIR
        yaml.dump(config, open(config_fpath, "w"), default_flow_style=False)
        data_dir = DEFAULT_DATA_DIR

    if not Path(data_dir).exists():
        raise OSError(f"Data directory specified in {data_dir} does not exist")

    return Path(data_dir)
