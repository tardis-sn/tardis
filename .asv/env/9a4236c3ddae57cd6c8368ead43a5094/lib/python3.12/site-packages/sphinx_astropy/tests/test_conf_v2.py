import pytest

pytest.importorskip("pydata_sphinx_theme")

from sphinx.cmd.build import build_main  # noqa: E402

BASIC_CONF = """
from sphinx_astropy.conf.v2 import *
suppress_warnings = ['app.add_directive', 'app.add_node', 'app.add_role']
"""

BASIC_INDEX = """
Title
=====

Just a test
"""


def generate_files(tmp_path):
    f1 = tmp_path / "conf.py"
    f1.write_text(BASIC_CONF)

    f2 = tmp_path / "index.rst"
    f2.write_text(BASIC_INDEX)


def test_conf(tmp_path):

    # Make sure the docs build with the v2 sphinx-astropy configuration

    generate_files(tmp_path)

    src_dir = str(tmp_path)
    html_dir = str(tmp_path / "html")

    status = build_main(argv=['-W', '-b', 'html', src_dir, html_dir])

    assert status == 0


def test_intersphinx_toggle(tmp_path, capsys):

    # Test the sphinx_astropy.ext.intersphinx_toggle extension

    generate_files(tmp_path)

    src_dir = str(tmp_path)
    html_dir = str(tmp_path / "html")

    status = build_main(argv=['-W', '-b', 'html', src_dir, html_dir, '-D',
                              'disable_intersphinx=1'])

    assert status == 0

    captured = capsys.readouterr()
    assert 'disabling intersphinx' in captured.out
    assert 'loading intersphinx' not in captured.out
