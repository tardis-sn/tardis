"""
    sphinxcontrib.apidoc.ext
    ~~~~~~~~~~~~~~~~~~~~~~~~

    A Sphinx extension for running 'sphinx-apidoc' on each build.

    :copyright: Copyright 2018 by Stephen Finucane <stephen@that.guru>
    :license: BSD, see LICENSE for details.
"""

from os import path

import sphinx
from sphinx.util import logging

try:
    from sphinx.ext import apidoc  # Sphinx >= 1.7
    _ignore_first_arg = False
except ImportError:
    from sphinx import apidoc  # Sphinx < 1.7
    _ignore_first_arg = True

if False:
    # For type annotation
    from sphinx.application import Sphinx  # noqa

logger = logging.getLogger(__name__)


def builder_inited(app):
    # type: (Sphinx) -> None
    module_dir = app.config.apidoc_module_dir
    output_dir = path.join(app.srcdir, app.config.apidoc_output_dir)
    excludes = app.config.apidoc_excluded_paths
    separate_modules = app.config.apidoc_separate_modules
    toc_file = app.config.apidoc_toc_file
    module_first = app.config.apidoc_module_first
    extra_args = app.config.apidoc_extra_args

    if toc_file and sphinx.version_info < (1, 8, 0):
        logger.warning("'apidoc_toc_file' is only supported by Sphinx "
                       "1.8+; skipping API doc generation")
        return

    if not module_dir:
        logger.warning("No 'apidoc_module_dir' specified; skipping API doc "
                       "generation")
        return

    # if the path is relative, make it relative to the 'conf.py' directory
    if not path.isabs(module_dir):
        module_dir = path.abspath(path.join(app.srcdir, module_dir))

    if not path.exists(module_dir):
        logger.warning(
            "The path defined in 'apidoc_module_dir' does not "
            "exist; skipping API doc generation; %s", module_dir)
        return

    # refactor this module so that we can call 'recurse_tree' like a sane
    # person - at present there is way too much passing around of the
    # 'optparse.Value' instance returned by 'optparse.parse_args'
    def cmd_opts():
        if _ignore_first_arg:
            yield 'sphinxcontrib-apidoc'

        yield '--force'

        if separate_modules:
            yield '--separate'

        if isinstance(toc_file, bool) and toc_file is False:
            yield '--no-toc'
        elif toc_file:
            yield '--tocfile'
            yield toc_file

        if module_first:
            yield '--module-first'

        yield '--output-dir'
        yield output_dir

        for arg in extra_args:
            yield arg

        yield module_dir

        for exc in excludes:
            yield path.abspath(path.join(module_dir, exc))

    apidoc.main(list(cmd_opts()))
