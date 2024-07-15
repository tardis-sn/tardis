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


am_replacer_str = """
This comes before

.. automodapi:: sphinx_automodapi.tests.example_module.mixed
{options}

This comes after
"""

am_replacer_basic_expected = """
This comes before


sphinx_automodapi.tests.example_module.mixed Module
---------------------------------------------------

.. automodule:: sphinx_automodapi.tests.example_module.mixed

Functions
^^^^^^^^^

.. automodsumm:: sphinx_automodapi.tests.example_module.mixed
    :functions-only:
    :toctree: api

Classes
^^^^^^^

.. automodsumm:: sphinx_automodapi.tests.example_module.mixed
    :classes-only:
    :toctree: api

Class Inheritance Diagram
^^^^^^^^^^^^^^^^^^^^^^^^^

.. automod-diagram:: sphinx_automodapi.tests.example_module.mixed
    :private-bases:
    :parts: 1

This comes after
"""


def test_am_replacer_basic(tmpdir):
    """
    Tests replacing an ".. automodapi::" with the automodapi no-option
    template
    """

    with open(tmpdir.join('index.rst').strpath, 'w') as f:
        f.write(am_replacer_str.format(options=''))

    run_sphinx_in_tmpdir(tmpdir)

    with open(tmpdir.join('index.rst.automodapi').strpath) as f:
        result = f.read()

    assert result == am_replacer_basic_expected


am_replacer_repr_str = u"""
This comes before with spéciàl çhars

.. automodapi:: sphinx_automodapi.tests.example_module.mixed
{options}

This comes after
"""


@pytest.mark.parametrize('writereprocessed', [False, True])
def test_am_replacer_writereprocessed(tmpdir, writereprocessed):
    """
    Tests the automodapi_writereprocessed option
    """

    with open(tmpdir.join('index.rst').strpath, 'w', encoding='utf-8') as f:
        f.write(am_replacer_repr_str.format(options=''))

    run_sphinx_in_tmpdir(tmpdir, additional_conf={'automodapi_writereprocessed': writereprocessed})

    assert tmpdir.join('index.rst.automodapi').isfile() is writereprocessed


am_replacer_noinh_expected = """
This comes before


sphinx_automodapi.tests.example_module.mixed Module
---------------------------------------------------

.. automodule:: sphinx_automodapi.tests.example_module.mixed

Functions
^^^^^^^^^

.. automodsumm:: sphinx_automodapi.tests.example_module.mixed
    :functions-only:
    :toctree: api

Classes
^^^^^^^

.. automodsumm:: sphinx_automodapi.tests.example_module.mixed
    :classes-only:
    :toctree: api


This comes after
""".format(empty='')


def test_am_replacer_noinh(tmpdir):
    """
    Tests replacing an ".. automodapi::" with no-inheritance-diagram
    option
    """

    ops = ['', ':no-inheritance-diagram:']
    ostr = '\n    '.join(ops)

    with open(tmpdir.join('index.rst').strpath, 'w') as f:
        f.write(am_replacer_str.format(options=ostr))

    run_sphinx_in_tmpdir(tmpdir)

    with open(tmpdir.join('index.rst.automodapi').strpath) as f:
        result = f.read()

    assert result == am_replacer_noinh_expected


am_replacer_titleandhdrs_expected = """
This comes before


sphinx_automodapi.tests.example_module.mixed Module
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

.. automodule:: sphinx_automodapi.tests.example_module.mixed

Functions
*********

.. automodsumm:: sphinx_automodapi.tests.example_module.mixed
    :functions-only:
    :toctree: api

Classes
*******

.. automodsumm:: sphinx_automodapi.tests.example_module.mixed
    :classes-only:
    :toctree: api

Class Inheritance Diagram
*************************

.. automod-diagram:: sphinx_automodapi.tests.example_module.mixed
    :private-bases:
    :parts: 1


This comes after
"""


def test_am_replacer_titleandhdrs(tmpdir):
    """
    Tests replacing an ".. automodapi::" entry with title-setting and header
    character options.
    """

    ops = ['', ':headings: &*']
    ostr = '\n    '.join(ops)

    with open(tmpdir.join('index.rst').strpath, 'w') as f:
        f.write(am_replacer_str.format(options=ostr))

    run_sphinx_in_tmpdir(tmpdir)

    with open(tmpdir.join('index.rst.automodapi').strpath) as f:
        result = f.read()

    assert result == am_replacer_titleandhdrs_expected


def test_am_replacer_titleandhdrs_invalid(tmpdir, capsys):
    """
    Tests replacing an ".. automodapi::" entry with title-setting and header
    character options.
    """

    ops = ['', ':headings: &']
    ostr = '\n    '.join(ops)

    with open(tmpdir.join('index.rst').strpath, 'w') as f:
        f.write(am_replacer_str.format(options=ostr))

    run_sphinx_in_tmpdir(tmpdir, expect_error=True)

    stdout, stderr = capsys.readouterr()
    assert "Not enough headings (got 1, need 2), using default -^" in stderr


am_replacer_nomain_str = """
This comes before

.. automodapi:: sphinx_automodapi.tests.example_module.functions
    :no-main-docstr:

This comes after
"""

am_replacer_nomain_expected = """
This comes before


sphinx_automodapi.tests.example_module.functions Module
-------------------------------------------------------



Functions
^^^^^^^^^

.. automodsumm:: sphinx_automodapi.tests.example_module.functions
    :functions-only:
    :toctree: api


This comes after
""".format(empty='')


def test_am_replacer_nomain(tmpdir):
    """
    Tests replacing an ".. automodapi::" with "no-main-docstring" .
    """

    with open(tmpdir.join('index.rst').strpath, 'w') as f:
        f.write(am_replacer_nomain_str.format(options=''))

    run_sphinx_in_tmpdir(tmpdir)

    with open(tmpdir.join('index.rst.automodapi').strpath) as f:
        result = f.read()

    assert result == am_replacer_nomain_expected


am_replacer_skip_str = """
This comes before

.. automodapi:: sphinx_automodapi.tests.example_module.functions
    :skip: add
    :skip: subtract

This comes after
"""

am_replacer_skip_expected = """
This comes before


sphinx_automodapi.tests.example_module.functions Module
-------------------------------------------------------

.. automodule:: sphinx_automodapi.tests.example_module.functions

Functions
^^^^^^^^^

.. automodsumm:: sphinx_automodapi.tests.example_module.functions
    :functions-only:
    :toctree: api
    :skip: add,subtract


This comes after
""".format(empty='')


def test_am_replacer_skip(tmpdir):
    """
    Tests using the ":skip: option in an ".. automodapi::" .
    """

    with open(tmpdir.join('index.rst').strpath, 'w') as f:
        f.write(am_replacer_skip_str.format(options=''))

    run_sphinx_in_tmpdir(tmpdir)

    with open(tmpdir.join('index.rst.automodapi').strpath) as f:
        result = f.read()

    assert result == am_replacer_skip_expected


am_replacer_skip_stdlib_str = """
This comes before

.. automodapi:: sphinx_automodapi.tests.example_module.stdlib
    :skip: time
    :skip: Path

This comes after
"""


am_replacer_skip_stdlib_expected = """
This comes before


sphinx_automodapi.tests.example_module.stdlib Module
----------------------------------------------------

.. automodule:: sphinx_automodapi.tests.example_module.stdlib

Functions
^^^^^^^^^

.. automodsumm:: sphinx_automodapi.tests.example_module.stdlib
    :functions-only:
    :toctree: api
    :skip: time,Path


This comes after
""".format(empty='')


def test_am_replacer_skip_stdlib(tmpdir):
    """
    Tests using the ":skip:" option in an ".. automodapi::"
    that skips objects imported from the standard library.
    This is a regression test for #141
    """

    with open(tmpdir.join('index.rst').strpath, 'w') as f:
        f.write(am_replacer_skip_stdlib_str.format(options=''))

    run_sphinx_in_tmpdir(tmpdir)

    with open(tmpdir.join('index.rst.automodapi').strpath) as f:
        result = f.read()

    assert result == am_replacer_skip_stdlib_expected


am_replacer_include_stdlib_str = """
This comes before

.. automodapi:: sphinx_automodapi.tests.example_module.stdlib
    :include: add
    :allowed-package-names: pathlib, datetime, sphinx_automodapi

This comes after
"""

am_replacer_include_stdlib_expected = """
This comes before


sphinx_automodapi.tests.example_module.stdlib Module
----------------------------------------------------

.. automodule:: sphinx_automodapi.tests.example_module.stdlib

Functions
^^^^^^^^^

.. automodsumm:: sphinx_automodapi.tests.example_module.stdlib
    :functions-only:
    :toctree: api
    :skip: Path,time
    :allowed-package-names: pathlib,datetime,sphinx_automodapi


This comes after
""".format(empty='')


def test_am_replacer_include_stdlib(tmpdir):
    """
    Tests using the ":include: option in an ".. automodapi::"
    in the presence of objects imported from the standard library.
    """

    with open(tmpdir.join('index.rst').strpath, 'w') as f:
        f.write(am_replacer_include_stdlib_str.format(options=''))

    run_sphinx_in_tmpdir(tmpdir)

    with open(tmpdir.join('index.rst.automodapi').strpath) as f:
        result = f.read()

    assert result == am_replacer_include_stdlib_expected


am_replacer_include_str = """
This comes before

.. automodapi:: sphinx_automodapi.tests.example_module.functions
    :include: add
    :include: subtract

This comes after
"""

am_replacer_include_expected = """
This comes before


sphinx_automodapi.tests.example_module.functions Module
-------------------------------------------------------

.. automodule:: sphinx_automodapi.tests.example_module.functions

Functions
^^^^^^^^^

.. automodsumm:: sphinx_automodapi.tests.example_module.functions
    :functions-only:
    :toctree: api
    :skip: multiply


This comes after
""".format(empty='')


def test_am_replacer_include(tmpdir):
    """
    Tests using the ":include: option in an ".. automodapi::" .
    """

    with open(tmpdir.join('index.rst').strpath, 'w') as f:
        f.write(am_replacer_include_str.format(options=''))

    run_sphinx_in_tmpdir(tmpdir)

    with open(tmpdir.join('index.rst.automodapi').strpath) as f:
        result = f.read()

    assert result == am_replacer_include_expected


am_replacer_invalidop_str = """
This comes before

.. automodapi:: sphinx_automodapi.tests.example_module.functions
    :invalid-option:

This comes after
"""


def test_am_replacer_invalidop(tmpdir, capsys):
    """
    Tests that a sphinx warning is produced with an invalid option.
    """

    with open(tmpdir.join('index.rst').strpath, 'w') as f:
        f.write(am_replacer_invalidop_str.format(options=''))

    run_sphinx_in_tmpdir(tmpdir, expect_error=True)

    stdout, stderr = capsys.readouterr()
    assert "Found additional options invalid-option in automodapi." in stderr


am_replacer_cython_str = """
This comes before

.. automodapi:: apyhtest_eva.unit02
{options}

This comes after
"""

am_replacer_cython_expected = """
This comes before


apyhtest_eva.unit02 Module
--------------------------

.. automodule:: apyhtest_eva.unit02

Functions
^^^^^^^^^

.. automodsumm:: apyhtest_eva.unit02
    :functions-only:
    :toctree: api

This comes after
""".format(empty='')


def test_am_replacer_cython(tmpdir, cython_testpackage):  # noqa
    """
    Tests replacing an ".. automodapi::" for a Cython module.
    """

    with open(tmpdir.join('index.rst').strpath, 'w') as f:
        f.write(am_replacer_cython_str.format(options=''))

    run_sphinx_in_tmpdir(tmpdir)

    with open(tmpdir.join('index.rst.automodapi').strpath) as f:
        result = f.read()

    assert result == am_replacer_cython_expected
