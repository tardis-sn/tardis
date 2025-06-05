# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("tardis")
except PackageNotFoundError:
    # Fallback for when setuptools_scm can't determine version
    __version__ = '2024.8.14.0.dev41+gf3a644caf5.d20240815'

__all__ = ['__version__', 'test']

# ----------------------------------------------------------------------------

import sys
import warnings

# ----------------------------------------------------------------------------

if ("astropy.units" in sys.modules) or ("astropy.constants" in sys.modules):
    warnings.warn(
        "Astropy is already imported externally. Astropy should be imported"
        " after TARDIS."
    )
else:
    from astropy import astronomical_constants, physical_constants

    physical_constants.set("codata2014")
    astronomical_constants.set("iau2012")

# ----------------------------------------------------------------------------

from tardis.base import run_tardis
from tardis.io.util import yaml_load_file as yaml_load
