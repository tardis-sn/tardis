.. _new-testing_and_validation:

Testing And Regression Data
===========================

Guidance for tests and regression data.

Tests
-----

.. _new-how-to-guide-run-unit-tests:

How-To Guide: Run Unit Tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TARDIS uses unit tests, regression-backed tests, and integration-style tests.
Run tests from the active ``tardis`` environment so imports, optional
dependencies, and regression-data tooling match the development checkout.

Run the unit tests from the repository root:

.. code-block:: shell

   pytest tardis


Run one package's tests:

.. code-block:: shell

   pytest tardis/io/model/readers/tests


Run one test file:

.. code-block:: shell

   pytest tardis/tests/test_util.py


.. _new-explanation-good-and-bad-test-cases:

Explanation: Good and Bad Test Cases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A good test case is narrow, deterministic, and tied to a behavior that matters.
It should fail for the bug or missing feature it describes, and it should not
depend on unrelated runtime state.

Good test case patterns:

- Tests one behavior at a time.
- Uses the smallest realistic fixture or input.
- Has a clear assertion.
- Uses ``numpy.testing`` helpers for numerical arrays.
- Uses ``pandas.testing`` helpers such as ``assert_frame_equal`` and
  ``assert_series_equal`` for pandas DataFrame and Series comparisons.
- Names the expected behavior in the test name.

Example:

.. code-block:: python

   import numpy.testing as npt


   def test_ascii_density_reader_converts_velocity_units():
       time_model, velocity, _ = read_simple_ascii_density(...)

       npt.assert_allclose(velocity[3].value, 1.3e4 * 1e5)
       npt.assert_almost_equal(time_model.to(u.day).value, 1.0)


A bad test case is broad, fragile, or unclear. It may pass while the behavior is
wrong, or fail because of unrelated changes.

Poor test case patterns:

- Tests many unrelated features at once.
- Uses a full simulation when a small parser/unit test would work.
- Asserts only that code ran without checking the result.
- Uses exact equality for floating-point arrays.
- Depends on local paths, network access, or test ordering.

Poor example:

.. code-block:: python

   def test_model():
       result = run_everything()
       assert result is not None


.. _new-reference-testing-command-reference:

Reference: Testing Command Reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run unit tests:

.. code-block:: shell

   pytest tardis

Run tests with coverage report generation:

.. code-block:: shell

   pytest tardis \
       --tardis-regression-data=/path/to/tardis-regression-data/ \
       --cov=tardis \
       --cov-report=xml \
       --cov-report=html


Clone regression data:

.. code-block:: shell

   GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/tardis-sn/tardis-regression-data.git
   cd tardis-regression-data
   git lfs fetch
   git lfs checkout


Run regression tests:

.. code-block:: shell

   pytest tardis --tardis-regression-data=/path/to/tardis-regression-data/


Run regression tests for a path:

.. code-block:: shell

   pytest tardis/path/to/test_file_or_directory \
       --tardis-regression-data=/path/to/tardis-regression-data/


Generate reference data:

.. code-block:: shell

   pytest tardis --tardis-regression-data=/path/to/tardis-regression-data --generate-reference


.. _new-explanation-testing-strategy:

Explanation: Testing Strategy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TARDIS uses a mix of unit tests, regression-backed tests, and integration-style
tests. Unit tests check individual functions and run quickly enough to provide
immediate feedback after changes, while regression-backed and integration-style
tests protect behavior that depends on larger simulations and stored reference
data.

The tests use ``pytest`` and are based on the ``astropy-setup-helpers`` package.

Unit tests are run after suggested changes to TARDIS to maintain code quality
and prevent regressions.

Regression Data
---------------

.. _new-explanation-regression-data:

Explanation: Regression Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Regression data tests compare TARDIS output against stored reference data. The
data lives in the ``tardis-regression-data`` repository and may include formats
such as NPY, HDF5, and CSV.

The project uses Git LFS for regression data. The regression data repository
should be cloned outside the main ``tardis`` repository.

The original documentation notes that TARDIS migrated from ``tardis-refdata`` to
``tardis-regression-data``.

.. _new-reference-regression-data-push-error:

Reference: Regression Data Push Error
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If pushing regression data fails with an error similar to:

.. code-block:: text

   remote: error: GH008: Your push referenced at least 1 unknown Git LFS object:
   remote:     <file-hash>
   remote: Try to push them with 'git lfs push --all'.


Check the LFS endpoint:

.. code-block:: shell

   git lfs env


Expected endpoint:

.. code-block:: text

   Endpoint=https://tardis:tardis-2025-lfs@registry.moria.egr.msu.edu/repository/tardis-lfs/info/lfs (auth=basic)


If the endpoint is different, check ``.git/config`` and ``.lfsconfig``. Git
prioritizes ``.git/config`` over ``.lfsconfig``, so remove duplicates and make sure
``.lfsconfig`` matches the TARDIS Regression Data repository. See the
``git-lfs-config`` documentation:

https://github.com/git-lfs/git-lfs/blob/main/docs/man/git-lfs-config.adoc

.. _new-how-to-guide-run-regression-data-tests:

How-To Guide: Run Regression Data Tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Advanced tests require TARDIS regression data from:

https://github.com/tardis-sn/tardis-regression-data

.. raw:: html

   <span style="color:red">Deleted: removed repeated migration and repository format explanation already covered in the regression-data explanation.</span>


Clone the regression data outside the main ``tardis`` repository:

.. code-block:: shell

   GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/tardis-sn/tardis-regression-data.git
   cd tardis-regression-data
   git lfs fetch
   git lfs checkout


``GIT_LFS_SKIP_SMUDGE=1`` defers downloading large files until the explicit
``git lfs fetch`` command. This makes download progress clearer and shows which
files are being retrieved.

Run regression tests:

.. code-block:: shell

   pytest tardis --tardis-regression-data=/path/to/tardis-regression-data/


Run tests for a specific file or directory:

.. code-block:: shell

   pytest tardis/path/to/test_file_or_directory \
       --tardis-regression-data=/path/to/tardis-regression-data/


For example, to run model-reader tests with local regression data:

.. code-block:: shell

   pytest tardis/io/model/readers/tests --tardis-regression-data=/path/to/tardis-regression-data/


Some cases require updating the regression data. See
:ref:`Update Regression Data <new-how-to-guide-update-regression-data>`.

The tests workflow runs on pull requests and push events. To prevent leaking LFS
quota, tests are disabled on forks. If you need to run tests on a fork, run the
tests workflow on the ``master`` branch first. The LFS cache generated on ``master``
should be available in child branches. Check whether the cache was generated in
the ``Setup LFS`` step of the workflow run. The cache can also be found under the
Management section of the Actions tab.

.. _new-how-to-guide-update-regression-data:

How-To Guide: Update Regression Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Regression data tests run only when ``pytest`` is called with the
``--tardis-regression-data`` flag. These tests compare TARDIS output, mostly
arrays, against information stored in regression data files.

If regression data tests fail in a pull request, possible causes include:

A. There is a problem in your code.
B. Your code is correct, but the regression data is outdated.
C. The pipeline is broken.

If you suspect the regression data is outdated:

1. Activate the ``tardis`` environment.
2. Fork and clone the ``tardis-regression-data`` repository.
3. Follow any instructions inside your local copy of ``tardis-regression-data``.
4. Go to your local ``tardis`` repository.
5. Make sure you are on the branch from which you want to generate new
   regression data.
6. Generate new regression data:

   .. code-block:: shell

      pytest tardis --tardis-regression-data=/path/to/tardis-regression-data --generate-reference


7. Check the results and make sure everything is correct.
8. Make a new branch in ``tardis-regression-data``.
9. Push the new regression data.
10. Open a pull request.

If issues arise, tag a TARDIS team member responsible for CI/CD:

https://tardis-sn.github.io/people/collaboration/

The pull-request regression comparison workflow uses the same idea in
``.github/workflows/compare-regdata.yml``:

.. code-block:: shell

   pytest tardis ${PYTEST_FLAGS} --generate-reference -m "not continuum"

.. _new-test-coverage:

Test Coverage
-------------

.. _new-how-to-guide-generate-and-view-test-coverage:

How-To Guide: Generate and View Test Coverage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TARDIS already uses coverage in CI. The test workflows pass coverage flags such
as ``--cov=tardis``, ``--cov-report=xml``, and ``--cov-report=html``, then combine
coverage artifacts with ``coverage combine``, ``coverage xml``, and ``coverage html``.
The older ``tox.ini`` coverage environment also writes ``coverage.xml``.

Run coverage from the active ``tardis`` environment. For TARDIS, coverage is only
meaningful when tests that depend on regression data are given the regression
data path; otherwise important regression-backed code paths may be skipped and
the coverage number is not a useful metric.

To generate local coverage for a targeted area:

.. code-block:: shell

   pytest tardis/io/model/readers/tests \
       --tardis-regression-data=/path/to/tardis-regression-data/ \
       --cov=tardis \
       --cov-report=xml \
       --cov-report=html


Then inspect the generated report:

.. code-block:: shell

   python -m http.server --directory htmlcov 8000


Open ``http://localhost:8000`` in a browser to inspect line-by-line coverage.

VS Code's built-in Testing view can also discover pytest tests, run or debug
tests, and run tests with coverage when the selected interpreter has the Python
testing dependencies installed. In a VS Code project, add the regression-data
argument to ``.vscode/settings.json`` so the Testing view runs TARDIS tests with
the same regression data used by the command line:

.. code-block:: json

   {
     "python.testing.pytestEnabled": true,
     "python.testing.pytestArgs": [
       "tardis",
       "--tardis-regression-data=[PATH_TO_REGRESSION_DATA]"
     ]
   }

After that, use the Testing view to discover tests and choose the coverage
action. VS Code can show coverage in the Test Coverage view and in editor
annotations when coverage data is available.

For VS Code, a popular option is the Coverage Gutters extension
(``ryanluker.vscode-coverage-gutters``). It displays coverage from XML or LCOV
files in the editor gutter.

Install it from VS Code Quick Open:

.. code-block:: text

   ext install ryanluker.vscode-coverage-gutters


Basic workflow:

1. Generate ``coverage.xml`` with pytest coverage.
2. Open the TARDIS repository in VS Code.
3. Run ``Coverage Gutters: Watch`` from the command palette.
4. Open a Python file. Covered, uncovered, and partially covered lines appear
   in the gutter.
5. Re-run tests and refresh coverage when the report changes.

See:

- VS Code testing documentation:
  https://code.visualstudio.com/docs/debugtest/testing
- Coverage Gutters marketplace page:
  https://marketplace.visualstudio.com/items?itemName=ryanluker.vscode-coverage-gutters
- Coverage.py command reference:
  https://coverage.readthedocs.io/en/latest/cmd.html
