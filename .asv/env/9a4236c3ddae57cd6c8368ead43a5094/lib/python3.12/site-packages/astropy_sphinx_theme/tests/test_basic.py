from distutils.version import LooseVersion

import pytest
from sphinx import __version__

SPHINX_LT_17 = LooseVersion(__version__) < LooseVersion('1.7')

if SPHINX_LT_17:
    from sphinx import build_main
else:
    from sphinx.cmd.build import build_main

THEMES = ['bootstrap-astropy']


BASIC_CONF = """
source_suffix = '.rst'
master_doc = 'index'
html_theme = '{theme}'
"""

BASIC_INDEX = """
Title
=====

Just a test
"""


@pytest.mark.parametrize('theme', THEMES)
def test_basic(tmpdir, theme):
    # Just make sure the docs build with the specified theme
    # (to make sure e.g. that no templates are missing)

    with open(tmpdir.join('conf.py').strpath, 'w') as f:
        f.write(BASIC_CONF.format(theme=theme))

    with open(tmpdir.join('index.rst').strpath, 'w') as f:
        f.write(BASIC_INDEX.format(theme=theme))

    src_dir = tmpdir.strpath
    html_dir = tmpdir.mkdir('html').strpath

    if SPHINX_LT_17:
        status = build_main(argv=['sphinx-build', '-W', '-b', 'html', src_dir, html_dir])
    else:
        # As of Sphinx 1.7, the first argument is now no longer ignored
        status = build_main(argv=['-W', '-b', 'html', src_dir, html_dir])

    assert status == 0
