Continuous Integration
----------------------

.. _explanation-continuous-integration:

Explanation: Continuous Integration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TARDIS uses continuous integration. When a change is proposed by pull request or
merged into ``master``, a service clones the repository, checks out the current
commit, and runs TARDIS tests. This helps detect bugs immediately.

TARDIS currently uses GitHub Actions for pipelines. A workflow is a YAML
configuration file with sections such as variables, jobs, and steps. Workflows
run commands when triggered by events, such as pushes or pull requests.

Making changes to an existing pipeline is done through a pull request. Creating
a new workflow requires adding a YAML file under ``.github/workflows``.

.. _how-to-guide-configure-a-github-actions-workflow:

How-To Guide: Configure a GitHub Actions Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To create a workflow, add a YAML file under ``.github/workflows``.


Use common TARDIS setup actions and settings:

1. Use the ``setup_lfs`` action and ``lfs-cache`` workflow when regression or atomic
   data is needed.
2. Use the ``lfs-cache`` workflow to cache regression data and atomic data and to
   check whether the cache is available.
3. Use the ``setup_env`` action to configure variables and settings for the
   pipeline.
4. Set the shell in the YAML configuration:

   .. code-block:: yaml

      defaults:
        run:
          shell: bash -l {0}

For workflows that inspect pull-request code:

1. Use ``pull_request`` and do not expose repository secrets to the job.
2. If the workflow must publish results or post authenticated comments, upload
   the results as artifacts and handle those operations in a separate trusted
   ``workflow_run`` workflow.
3. Use ``pull_request_target`` only for checkout-free tasks such as labeling or
   contributor notifications.

The ``docs``, ``codestyle``, ``mailmap``, and ORCID workflows use this split
between unprivileged checks and trusted publishers. The utility workflow keeps
``pull_request_target`` only for its checkout-free label and
first-contributor actions.


The ``docs`` workflow in ``.github/workflows/build-docs.yml`` uses the shared LFS
cache workflow, restores sparse atomic data, sets up the environment, and builds
the docs:

.. code-block:: yaml

   jobs:
     test-cache:
       uses: ./.github/workflows/lfs-cache.yml
       with:
         atom-data-sparse: true
         regression-data-repo: tardis-sn/tardis-regression-data

     build-docs:
       steps:
         - name: Setup LFS
           uses: ./.github/actions/setup_lfs
           with:
             atom-data-sparse: true
         - name: Setup environment
           uses: tardis-sn/tardis-actions/setup-env@main
         - name: Build documentation
           run: cd docs/ && make html NCORES=auto

.. _reference-continuous-integration-reference:

Reference: Continuous Integration Reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Cache Keys
^^^^^^^^^^

Regression data cache key format:

.. code-block:: text

   tardis-regression-<data-type>-<hash>-v1


Examples:

- ``tardis-regression-atom-data-sparse-<hash>-v1``: atomic data cache.
- ``tardis-regression-full-data-<hash>-v1``: full TARDIS regression data cache.

Used in the ``setup_lfs`` action.

Environment cache key format:

.. code-block:: text

   tardis-conda-env-<os-label>-<hash>-v1


Examples:

- ``tardis-conda-env-linux-<hash>-v1``: Linux conda environment.
- ``tardis-conda-env-macos-<hash>-v1``: macOS conda environment.

Used in the ``setup_env`` action.

The ``-v1`` suffix allows future cache invalidation. The ``lfs-cache`` workflow
fails if the cache is unavailable and does not pull LFS data by default. If the
``allow_lfs_pull`` label is added to a pull request, the workflow pulls LFS data.
Use this label sparingly and only with caution.

Common Pipeline Steps
^^^^^^^^^^^^^^^^^^^^^

- Use ``setup_lfs`` and the ``lfs-cache`` workflow when regression or atomic data is
  required.
- Use ``setup_env`` to configure environment variables and settings.
- Use ``bash -l {0}`` as the run shell in workflow YAML.

Documentation Build Pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The documentation build pipeline builds and deploys the TARDIS documentation
website.

Documentation Preview Pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The documentation preview pipeline does not run on the main repository, only on
forks. It supports pull request documentation previews.

Testing Pipeline
^^^^^^^^^^^^^^^^

The testing pipeline runs four concurrent jobs:

- Ubuntu tests without the continuum marker.
- Ubuntu tests with the continuum marker.
- macOS tests without the continuum marker.
- macOS tests with the continuum marker.

This is two platforms multiplied by two test types. The pipeline includes
environment installation, regression data configuration, and coverage report
upload after tests finish.

GPU tests can be triggered by applying the ``full-tests`` label to a pull request.

Release Pipeline
^^^^^^^^^^^^^^^^

TARDIS publishes a new release every Sunday at 00:00 UTC.

Pre-Release
^^^^^^^^^^^

The pre-release action clones ``tardis-sn/tardis_zenodo``, runs a notebook to
generate ``.zenodo.json``, and pushes the file to the root of the TARDIS
repository. The file creates a new TARDIS version on Zenodo with all committers
as authors. A pull request is created and automatically merged if required
checks pass.

Zenodo job:

1. Check out ``tardis-sn/tardis_zenodo``.
2. Wait for the Zenodo webhook to be available with a 3 minute sleep.
3. Set up the Python environment stored in ``tardis-sn/tardis_zenodo``.
4. Store the Zenodo API secret key in an environment variable.
5. Run the notebook to generate ``.zenodo.json``.
6. Re-run if there are errors and ignore errors.
7. Upload ``.zenodo.json`` as an artifact.

Pip tests job:

Runs the TARDIS test suite after installing TARDIS from ``master`` with pip using
Git, instead of an editable install. This catches cases where a new module lacks
an ``__init__.py`` file or data files are missing from
``[tool.setuptools.package-data]`` in ``pyproject.toml``.

Pre-release pull request job:

1. Relies on Zenodo and pip test steps completing.
2. Checks out the TARDIS repository.
3. Downloads artifacts from previous steps.
4. Checks for ``.zenodo.json`` and uses it if generated.
5. Gets the current date.
6. Creates a bot pull request on the ``tardis-bot`` fork using branch
   ``pre-release-<date>`` with the new ``.zenodo.json``.
7. Waits 1 minute for the pull request to be created.
8. Automatically approves the pull request using tokens from infrastructure and
   core coordinator members.
9. Enables auto-merge.

Release
^^^^^^^

The release job creates a GitHub release after the pre-release pull request is
merged.

1. Check out TARDIS with fetch depth 0.
2. Set up Python.
3. Install ``setuptools_scm`` and ``git-cliff``.
4. Get the current TARDIS version using ``setuptools_scm`` via a helper script.
5. Get the next TARDIS version using ``setuptools_scm``.
6. Create a GitHub release using the new version as the tag.
7. Wait 2 minutes for Zenodo to update the new TARDIS release.
8. Fetch the new DOI from Zenodo using the Zenodo API and create a badge.
9. Generate the changelog with ``git-cliff``.
10. Update the release description with the changelog and Zenodo badge.
11. Include environment lock files in release assets.

Post-Release
^^^^^^^^^^^^

The post-release action updates the changelog, citation, and credits in the main
repository.

Changelog job:

1. Check out TARDIS with fetch depth 0.
2. Get the current release tag.
3. Generate a changelog with ``git-cliff``.
4. Upload ``CHANGELOG.md`` as an artifact.

Citation job:

1. Check out TARDIS.
2. Wait 3 minutes for the Zenodo webhook.
3. Set up Python.
4. Install ``doi2cff``.
5. Convert the latest TARDIS release DOI to ``CITATION.cff``.
6. Try 10 times with a 60 second sleep between attempts.
7. Upload ``CITATION.cff`` as an artifact.

Credits job:

1. Check out TARDIS.
2. Wait 3 minutes for the Zenodo webhook.
3. Set up Python.
4. Install ``requests``.
5. Run a helper script to update ``README.rst`` and ``docs/resources/credits.rst``.
6. Upload ``README.rst`` and ``credits.rst`` as artifacts.
7. Dispatch updates to the TARDIS website.

Post-release pull request job:

1. Check out TARDIS.
2. Download artifacts from previous steps.
3. Copy ``CHANGELOG.md``, ``CITATION.cff``, ``README.rst``, and ``credits.rst`` into
   the repository.
4. Get the current date.
5. Create a pull request.
6. Wait 30 seconds for the pull request to be created.
7. Automatically approve the pull request using tokens from infrastructure and
   core coordinator members.
8. Enable auto-merge.

TARDIS Carsus Compatibility Check
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The TARDIS Carsus Compatibility Check, or Bridge, compares reference data
generated with different Carsus versions. It has two jobs:

- ``carsus-build``: generates an atomic file with the latest Carsus.
- ``tardis-build``: generates new reference data with the atomic file.

The two reference data files are compared with this notebook:

https://github.com/tardis-sn/tardis-refdata/blob/master/notebooks/ref_data_compare_from_paths.ipynb

The workflow has a ``workflow_dispatch`` event for manual runs and is also
triggered weekly by the save atomic files workflow.

Save Atomic Files Workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^

The Save Atomic Files workflow runs weekly and can be triggered manually. It
runs the Bridge and sends an artifact to Moria containing the generated atomic
data file and comparison notebook. A separate job indicates whether the Bridge
failed.

Regression Data Comparison Workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Regression Data Comparison workflow compares regression data between the
current branch and the base branch on pull requests. It runs only on pull
requests, not on ``master``.

The workflow generates regression data for the latest commit on the pull request
and compares it with ``master`` using the comparison notebook. The notebook is
uploaded as an artifact and pushed to the ``reg-data-comp`` repository for
previews in the bot comment.

The comparison is split into two workflows for security. The
``compare-regdata`` workflow runs pull-request code without repository secrets
and uploads the comparison artifacts. The ``publish-regdata-comparison``
workflow runs after a successful comparison, checks the pull-request label, and
publishes those artifacts and the bot comment with its trusted credentials.

The workflow exports images from the comparison notebook and embeds them in the
bot comment. Unless there are key changes to HDF files in the regression data,
the bot shows two images: spectrum change and relative changes in keys. If
there are key changes, the bot shows a third image visualizing key changes.

LFS-Cache Workflow
^^^^^^^^^^^^^^^^^^

The ``LFS-cache`` workflow caches regression data and atomic data. It can be
triggered manually or when there is a push to the main branch of the regression
data repository. It performs LFS pulls when necessary and caches objects. The
``setup-lfs`` action restores cached objects. Both fail if the cache is
unavailable.
