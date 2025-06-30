# Licensed under a 3-clause BSD style license - see LICENSE.rst
import sys

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from importlib.metadata import version as ilversion
from packaging.version import Version as pversion

__version__ = ilversion("tardis")
last_release = pversion(__version__).base_version
__all__ = ['__version__', 'run_tardis', 'last_release']

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
