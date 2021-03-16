# Licensed under a 3-clause BSD style license - see LICENSE.rst
import warnings

import pyne.data

import tardis.util.custom_logger as custom_logger_settings

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *

# ----------------------------------------------------------------------------

from tardis.base import run_tardis
from tardis.io.util import yaml_load_config_file as yaml_load

warnings.filterwarnings("ignore", category=pyne.utils.QAWarning)

custom_logger_settings.init()
custom_logger_settings.save = False

# the default level is TARDIS INFO
# due to this, no WARNINGS are shown by default
custom_logger_settings.level = "TARDIS INFO"
custom_logger_settings.reset_logger()
