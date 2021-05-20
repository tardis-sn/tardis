# Licensed under a 3-clause BSD style license - see LICENSE.rst
import sys
import logging
import warnings

import pyne.data

from tardis.util.colored_logger import ColoredFormatter, formatter_message

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *

# ----------------------------------------------------------------------------

from tardis.base import run_tardis
from tardis.io.util import yaml_load_config_file as yaml_load

warnings.filterwarnings("ignore", category=pyne.utils.QAWarning)

FORMAT = "[$BOLD%(name)-20s$RESET][%(levelname)-18s]  %(message)s ($BOLD%(filename)s$RESET:%(lineno)d)"
COLOR_FORMAT = formatter_message(FORMAT, True)

logging.captureWarnings(True)
logger = logging.getLogger("tardis")
logger.setLevel(logging.INFO)

console_handler = logging.StreamHandler(sys.stdout)
console_formatter = ColoredFormatter(COLOR_FORMAT)
console_handler.setFormatter(console_formatter)

logger.addHandler(console_handler)
logging.getLogger("py.warnings").addHandler(console_handler)

# ----------------------------------------------------------------------------
# pyne holds Python 3.7 on macOS, but refdata is pickled with protocol 5 (3.8.3)

if sys.version_info < (3, 8, 3):
    import pickle5

    sys.modules["pickle"] = pickle5
