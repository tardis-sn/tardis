# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from copy import copy

import pytest
from docutils.parsers.rst import directives, roles

from . import cython_testpackage  # noqa
from .helpers import run_sphinx_in_tmpdir


def setup_function(func):
    # This can be replaced with the docutils_namespace context manager once
    # it is in a stable release of Sphinx
    func._directives = copy(directives._directives)
    func._roles = copy(roles._roles)


def teardown_function(func):
    directives._directives = func._directives
    roles._roles = func._roles


# nosignatures

ADD_RST = """
:orphan:

add
===

.. currentmodule:: sphinx_automodapi.tests.example_module.mixed

.. autofunction:: add
""".strip()

MIXEDSPAM_RST = """
:orphan:

MixedSpam
=========

.. currentmodule:: sphinx_automodapi.tests.example_module.mixed

.. autoclass:: MixedSpam
   :show-inheritance:
""".strip()


def write_api_files_to_tmpdir(tmpdir):
    apidir = tmpdir.mkdir('api')
    with open(apidir.join('sphinx_automodapi.tests.example_module.mixed.add.rst').strpath, 'w') as f:
        f.write(ADD_RST)
    with open(apidir.join('sphinx_automodapi.tests.example_module.mixed.MixedSpam.rst').strpath, 'w') as f:
        f.write(MIXEDSPAM_RST)


ams_to_asmry_str = """
Before

.. automodsumm:: sphinx_automodapi.tests.example_module.mixed
{options}

And After
"""

ams_to_asmry_expected = """\
.. currentmodule:: sphinx_automodapi.tests.example_module.mixed

.. autosummary::

    add
    MixedSpam

"""


def test_ams_to_asmry(tmpdir):

    with open(tmpdir.join('index.rst').strpath, 'w') as f:
        f.write(ams_to_asmry_str.format(options=''))

    write_api_files_to_tmpdir(tmpdir)

    run_sphinx_in_tmpdir(tmpdir)

    with open(tmpdir.join('index.rst.automodsumm').strpath) as f:
        result = f.read()

    assert result == ams_to_asmry_expected


def test_too_many_options(tmpdir, capsys):

    ops = ['', ':classes-only:', ':functions-only:']
    ostr = '\n    '.join(ops)

    with open(tmpdir.join('index.rst').strpath, 'w') as f:
        f.write(ams_to_asmry_str.format(options=ostr))

    write_api_files_to_tmpdir(tmpdir)

    run_sphinx_in_tmpdir(tmpdir, expect_error=True)

    stdout, stderr = capsys.readouterr()
    assert ("[automodsumm] Defined more than one of functions-only, "
            "classes-only, and variables-only.  Skipping this directive." in stderr)


ORDEREDDICT_RST = """
:orphan:

OrderedDict
===========

.. currentmodule:: sphinx_automodapi.tests.example_module.noall

.. autoclass:: OrderedDict
   :show-inheritance:
""".strip()


@pytest.mark.parametrize('options,expect', [
    ('', ['add', 'MixedSpam']),
    (':allowed-package-names: sphinx_automodapi', ['add', 'MixedSpam']),
    (':allowed-package-names: collections', ['OrderedDict']),
    (':allowed-package-names: sphinx_automodapi,collections',
     ['add', 'MixedSpam', 'OrderedDict']),
])
def test_am_allowed_package_names(options, expect, tmpdir):
    """
    Test that allowed_package_names is interpreted correctly.
    """
    def mixed2noall(s):
        return s.replace('example_module.mixed', 'example_module.noall')

    am_str = ams_to_asmry_str
    with open(tmpdir.join('index.rst').strpath, 'w') as f:
        f.write(mixed2noall(am_str).format(options=('   '+options if options else '')))

    apidir = tmpdir.mkdir('api')
    with open(apidir.join('sphinx_automodapi.tests.example_module.noall.add.rst').strpath, 'w') as f:
        f.write(mixed2noall(ADD_RST))
    with open(apidir.join('sphinx_automodapi.tests.example_module.noall.MixedSpam.rst').strpath, 'w') as f:
        f.write(mixed2noall(MIXEDSPAM_RST))
    with open(apidir.join('sphinx_automodapi.tests.example_module.noall.OrderedDict.rst').strpath, 'w') as f:
        f.write(ORDEREDDICT_RST)

    run_sphinx_in_tmpdir(tmpdir)

    with open(tmpdir.join('index.rst.automodsumm').strpath) as f:
        result = f.read()

    for x in expect:
        assert '   '+x in result


PILOT_RST = """
:orphan:

pilot
=====

.. currentmodule:: apyhtest_eva.unit02

.. autofunction:: pilot
""".strip()

ams_cython_str = """
Before

.. automodsumm:: apyhtest_eva.unit02
    :functions-only:

And After
"""

ams_cython_expected = """\
.. currentmodule:: apyhtest_eva.unit02

.. autosummary::

    pilot

"""


def test_ams_cython(tmpdir, cython_testpackage):  # noqa

    with open(tmpdir.join('index.rst').strpath, 'w') as f:
        f.write(ams_cython_str)

    apidir = tmpdir.mkdir('api')
    with open(apidir.join('apyhtest_eva.unit02.pilot.rst').strpath, 'w') as f:
        f.write(PILOT_RST)

    run_sphinx_in_tmpdir(tmpdir)

    with open(tmpdir.join('index.rst.automodsumm').strpath) as f:
        result = f.read()

    assert result == ams_cython_expected
