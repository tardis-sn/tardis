"""
    sphinxcontrib.apidoc
    ~~~~~~~~~~~~~~~~~~~~

    A Sphinx extension for running 'sphinx-apidoc' on each build.

    :copyright: Copyright 2018 by Stephen Finucane <stephen@that.guru>
    :license: BSD, see LICENSE for details.
"""

import pbr.version

from sphinxcontrib.apidoc import ext

__version__ = pbr.version.VersionInfo('sphinxcontrib-apidoc').version_string()

if False:
    # For type annotation
    from typing import Any, Dict  # noqa
    from sphinx.application import Sphinx  # noqa


def setup(app):
    # type: (Sphinx) -> Dict[unicode, Any]
    app.setup_extension('sphinx.ext.autodoc')  # We need autodoc to function

    app.connect('builder-inited', ext.builder_inited)
    app.add_config_value('apidoc_module_dir', None, 'env', [str])
    app.add_config_value('apidoc_output_dir', 'api', 'env', [str])
    app.add_config_value('apidoc_excluded_paths', [], 'env', [[str]])
    app.add_config_value('apidoc_separate_modules', False, 'env', [bool])
    app.add_config_value('apidoc_toc_file', None, 'env', [str, bool])
    app.add_config_value('apidoc_module_first', False, 'env', [bool])
    app.add_config_value('apidoc_extra_args', [], 'env', [list])

    return {'version': __version__, 'parallel_read_safe': True}
