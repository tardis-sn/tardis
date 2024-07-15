# This directory contains the default Sphinx configuration for the core Astropy
# package and other packages that want to use the same configuration. Rather
# than store the configuration in a single conf.py file, we store it in v?.py
# files and import the latest one into sphinx_astropy.conf so that packages can
# choose to use either to do:
#
#     from sphinx_astropy.conf import *
#
# or:
#
#     from sphinx_astropy.conf.v1 import *
#
# with the latter being the option to use for stability. The idea is that
# we can still make small changes (mainly fixing bugs) to v1.py, but if we
# make any big changes in future, we can create a new version that packages
# can choose to opt-in to. To create a new default configuration, create a
# v2.py file (either starting from a copy of v1.py or starting from
# scratch), and change the import below to 'from .v2 import *'.

# TODO: Switch default to v2
from .v1 import *  # noqa: F401, F403
