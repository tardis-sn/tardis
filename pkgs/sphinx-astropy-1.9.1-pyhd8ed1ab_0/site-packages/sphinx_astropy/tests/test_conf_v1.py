from sphinx.cmd.build import build_main

BASIC_CONF = """
from sphinx_astropy.conf import *
suppress_warnings = ['app.add_directive', 'app.add_node', 'app.add_role']
"""

BASIC_INDEX = """
Title
=====

Just a test
"""


def generate_files(tmpdir):

    with open(tmpdir.join('conf.py').strpath, 'w') as f:
        f.write(BASIC_CONF)

    with open(tmpdir.join('index.rst').strpath, 'w') as f:
        f.write(BASIC_INDEX)


def test_conf(tmpdir):

    # Just make sure the docs build with the default sphinx-astropy configuration

    generate_files(tmpdir)

    src_dir = tmpdir.strpath
    html_dir = tmpdir.mkdir('html').strpath

    status = build_main(argv=['-W', '-b', 'html', src_dir, html_dir])

    assert status == 0


def test_intersphinx_toggle(tmpdir, capsys):

    # Test the sphinx_astropy.ext.intersphinx_toggle extension

    generate_files(tmpdir)

    src_dir = tmpdir.strpath
    html_dir = tmpdir.mkdir('html').strpath

    status = build_main(argv=['-W', '-b', 'html', src_dir, html_dir, '-D', 'disable_intersphinx=1'])

    assert status == 0

    captured = capsys.readouterr()
    assert 'disabling intersphinx' in captured.out
    assert 'loading intersphinx' not in captured.out
