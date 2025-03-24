# Licensed under a 3-clause BSD style license - see LICENSE.rst
import sys

__all__ = ['__version__', 'run_tardis', 'yaml_load']

try:
    from .version import version as __version__
except ImportError:
    __version__ = ''

if ("astropy.units" in sys.modules) or ("astropy.constants" in sys.modules):
    import warnings

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
