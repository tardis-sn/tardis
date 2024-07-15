# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
The purpose of this extension is to provide a configuration value that can be
used to disable intersphinx on the command-line without editing conf.py. To use,
you can build documentation with::

    sphinx-build ... -D disable_intersphinx=1

This is used e.g. by astropy-helpers when using the build_docs command.
"""

from sphinx.util.console import bold
from sphinx.util import logging


def disable_intersphinx(app, config=None):
    info = logging.getLogger(__name__).info

    if app.config.disable_intersphinx:
        info(bold('disabling intersphinx...'))
        app.config.intersphinx_mapping.clear()


def setup(app):
    app.connect('config-inited', disable_intersphinx)
    app.add_config_value('disable_intersphinx', 0, True)
    return {'parallel_read_safe': True,
            'parallel_write_safe': True}
