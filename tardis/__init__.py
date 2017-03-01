# Licensed under a 3-clause BSD style license - see LICENSE.rst
import logging
import color
# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

from tardis.base import run_tardis
from tardis.io.util import yaml_load_config_file as yaml_load

logging.captureWarnings(True)
logger = logging.getLogger('tardis')
logging.setLoggerClass(color.ColoredLogger)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(console_formatter)
logger.addHandler(console_handler)
logging.getLogger('py.warnings').addHandler(console_handler)
