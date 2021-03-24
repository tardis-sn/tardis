from tardis import __path__ as TARDIS_PATH
import os, shutil
import yaml

from astropy.config import get_config_dir

TARDIS_PATH = TARDIS_PATH[0]
DEFAULT_CONFIG_PATH = os.path.join(
    TARDIS_PATH, "data", "default_tardis_internal_config.yml"
)
DEFAULT_DATA_DIR = os.path.join(
    os.path.expanduser("~"), "Downloads", "tardis-data"
)
from tardis.util.custom_logger import logger


def get_internal_configuration():

    config_fpath = os.path.join(get_config_dir(), "tardis_internal_config.yml")
    if not os.path.exists(config_fpath):
        logger.warning(
            "Configuration File {0} does not exist - creating new one from default".format(
                config_fpath
            )
        )
        shutil.copy(DEFAULT_CONFIG_PATH, config_fpath)
    with open(config_fpath) as config_fh:
        return yaml.load(config_fh, Loader=yaml.CLoader)


def get_data_dir():

    config = get_internal_configuration()
    data_dir = config.get("data_dir", None)
    if data_dir is None:
        config_fpath = os.path.join(
            get_config_dir(), "tardis_internal_config.yml"
        )
        logger.critical(
            "\n{line_stars}\n\nTARDIS will download different kinds of data (e.g. atomic) to its data directory {default_data_dir}\n\n"
            "TARDIS DATA DIRECTORY not specified in {config_file}:\n\n"
            "ASSUMING DEFAULT DATA DIRECTORY {default_data_dir}\n "
            "YOU CAN CHANGE THIS AT ANY TIME IN {config_file} \n\n"
            "{line_stars} \n\n".format(
                line_stars="*" * 80,
                config_file=config_fpath,
                default_data_dir=DEFAULT_DATA_DIR,
            )
        )
        if not os.path.exists(DEFAULT_DATA_DIR):
            os.makedirs(DEFAULT_DATA_DIR)
        config["data_dir"] = DEFAULT_DATA_DIR
        yaml.dump(config, open(config_fpath, "w"), default_flow_style=False)
        data_dir = DEFAULT_DATA_DIR

    if not os.path.exists(data_dir):
        raise IOError(
            "Data directory specified in {0} does not exist".format(data_dir)
        )

    return data_dir
