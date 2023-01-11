# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *  # noqa

# ----------------------------------------------------------------------------

__all__ = []

# ----------------------------------------------------------------------------

import sys

# ----------------------------------------------------------------------------

from astropy import physical_constants, astronomical_constants

physical_constants.set("codata2014")
astronomical_constants.set("iau2012")

# ----------------------------------------------------------------------------

from tardis.base import run_tardis
from tardis.io.util import yaml_load_file as yaml_load
