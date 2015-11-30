**************
Reporting Bugs
**************

Despite our best efforts TARDIS contains bugs. We want to make sure that fixing
these is quick and painless for users and developers.

When reporting bugs please use the mailing list
(`tardis-sn-user <http://groups.google.com/forum/#!forum/tardis-sn-users>`_) or
create an issue for TARDIS on github (??? ;preferred).

If there are very long ouputs needed to debug the problem, please use a website like
`pastebin <http://pastebin.com>`_ that way it is easier for us to look through it.


One of the major problems that we see are mismatches in installed versions of
third-party libraries (that's where virtualenvs really help -
see :doc:`workflow/python_environment`). To make sure that we can test this properly
please tell us the installed version of third party packages. The easiest way is::


    $ pip freeze
    Cython==0.20.2
    Jinja2==2.7.2
    MarkupSafe==0.21
    PyYAML==3.11
    Pygments==1.6
    Sphinx==1.2.2
    astropy==1.0.dev10065
    astropy-helpers==0.4.2
    backports.ssl-match-hostname==3.4.0.2
    cov-core==1.11
    coverage==3.7.1
    docutils==0.11
    execnet==1.2.0
    h5py==2.4.0a0
    ipdb==0.8
    ipython==1.2.1
    latexcodec==0.3.2
    matplotlib==1.4.0
    mock==1.0.1
    nose==1.3.4
    numexpr==2.4
    numpy==1.9.0
    numpydoc==0.4
    oset==0.1.3
    pandas==0.14.1
    pep8==1.5.6
    py==1.4.20
    pybtex==0.17
    pybtex-docutils==0.2.0
    pyparsing==2.0.2
    pytest==2.5.2
    pytest-cache==1.0
    pytest-cov==1.6
    pytest-ipdb==0.1-prerelease
    pytest-pep8==1.0.6
    python-dateutil==2.2
    pytz==2014.7
    pyzmq==14.2.0
    scipy==0.13.3
    six==1.8.0
    sphinx-bootstrap-theme==0.4.0
    sphinxcontrib-bibtex==0.2.9
    sphinxcontrib-tikz==0.4.1
    tables==3.1.1
    tornado==3.2
    tox==1.7.1
    virtualenv==1.11.6
    wsgiref==0.1.2


and pasting the output in your error report, this will show us what packages are
installed and used.

The next thing that is useful to do if TARDIS is downloaded in a directory is to run
the test suites and find if any of the tests fail::

    python setup.py test
    running test
    running build
    running build_py
    copying tardis/./macro_atom.c -> build/lib.macosx-10.9-x86_64-2.7/tardis/.
    copying tardis/./montecarlo.c -> build/lib.macosx-10.9-x86_64-2.7/tardis/.
    running build_ext
    skipping 'tardis/montecarlo.c' Cython extension (up-to-date)
    skipping 'tardis/macro_atom.c' Cython extension (up-to-date)
    running build_scripts
    ============================= test session starts ==============================
    platform darwin -- Python 2.7.8 -- pytest-2.5.1

    Running tests in tardis /Users/wkerzend/scripts/python/tardis/docs.

    Platform: Darwin-13.4.0-x86_64-i386-64bit

    Executable: /Users/wkerzend/.virtualenvs/tardis-devel/bin/python

    Full Python Version:
    2.7.8 (default, Oct 15 2014, 22:04:42)
    [GCC 4.2.1 Compatible Apple LLVM 5.1 (clang-503.0.40)]

    encodings: sys: ascii, locale: UTF-8, filesystem: utf-8, unicode bits: 15
    byteorder: little
    float info: dig: 15, mant_dig: 15

    Numpy: 1.9.1
    Scipy: 0.13.3
    astropy: 1.0.dev10065
    yaml: 3.11
    cython: 0.20.2
    h5py: 2.4.0a0
    Matplotlib: 1.4.0
    ipython: 1.2.1

    plugins: cache, cov, ipdb, pep8
    collected 107 items

    tardis/io/tests/test_ascii_readers.py .......
    tardis/io/tests/test_config_reader.py .............................
    tardis/io/tests/test_config_validator.py ........................
    tardis/io/tests/test_configuration_namespace.py .........
    tardis/tests/test_atomic.py .....s
    tardis/tests/test_lte_plasma.py ssssssssss
    tardis/tests/test_plasma_nlte.py .....
    tardis/tests/test_plasma_simple.py .
    tardis/tests/test_tardis_full.py s
    tardis/tests/test_util.py ...............

    ==================== 95 passed, 12 skipped in 8.53 seconds =====================



This will hopefully help us to identify us the problem. We will continue to
update this document with other techniques.