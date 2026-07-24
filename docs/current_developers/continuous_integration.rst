Continuous Integration
----------------------

.. _explanation-continuous-integration:

Explanation: Continuous Integration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TARDIS uses continuous integration. When a change is proposed by pull request,
pushed to a configured branch, or reached by a scheduled workflow, a service
clones the repository, checks out the relevant commit, and runs the appropriate
checks. This helps detect bugs immediately.

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
2. If the workflow must publish results or post authenticated comments, follow
   :ref:`secure-pull-request-publishing`.
3. Use ``pull_request_target`` only for checkout-free tasks such as labeling or
   contributor notifications, or check out only an explicitly trusted branch.

.. _secure-pull-request-publishing:

Secure pull-request publishing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TARDIS uses two workflow runs when a pull request needs both untrusted code
execution and a write operation. This is a security boundary, not merely an
organizational split:

1. The *producer* uses ``pull_request``. It may check out, install, test, or
   build the pull-request code, but it has no repository secrets.
2. The *publisher* uses ``workflow_run``. GitHub loads this workflow from the
   default branch, where it can use a write-capable token. It downloads output
   from the producer's run ID and SHA-specific artifact, but never checks out or
   executes pull-request code.

Keep artifacts at this boundary data-only. A publisher may copy generated HTML,
read a text log, or include an image in a comment. It must not run a script,
action, executable, or package taken from the artifact. See GitHub's
`pull_request_target security guidance <https://docs.github.com/en/actions/reference/security/securely-using-pull_request_target>`_
and the security warning in the
`workflow_run documentation <https://docs.github.com/en/actions/using-workflows/events-that-trigger-workflows#workflow_run>`_.

GitHub sometimes supplies an empty ``workflow_run.pull_requests`` array for a
fork pull request. This cannot be changed by using a different token, so do not
read the PR number from that array. TARDIS publishers instead use the
`Find Pull Request action <https://github.com/marketplace/actions/find-pull-request>`_
in a fail-closed ``Get PR number`` step:

.. code-block:: yaml

   - name: Get PR number
     id: pr
     uses: juliangruber/find-pull-request-action@v1
     with:
       github-token: ${{ github.token }}
       branch: ${{ github.event.workflow_run.head_repository.owner.login }}:${{ github.event.workflow_run.head_branch }}
       state: open

The action searches the current repository for open pull requests whose head is
the source repository and branch from the completed producer run. The
``owner:branch`` form is important for fork pull requests; using only the branch
name could match a branch with the same name from another repository.

Publishers guard their artifact, deployment, and comment steps with all of the
following action outputs:

* ``steps.pr.outputs.number`` is not empty;
* ``steps.pr.outputs.matching-pr-count`` is ``1``; and
* ``steps.pr.outputs.head-sha`` equals ``env.PR_SHA``, which is
  ``github.event.workflow_run.head_sha``.

These checks preserve the previous fail-closed behavior. They prevent an old
producer run from publishing after a newer commit has been pushed, and prevent
publication when the branch selector is ambiguous. The action handles API
authentication and exposes the pull-request number and head SHA as step
outputs, so the publishers no longer need a shell script or the GitHub CLI for
PR discovery.

If the action's API request fails, for example because the token lacks
permission, the ``Get PR number`` step fails instead of silently skipping
publication. The publisher job's top-level condition also ensures this step
runs only when the producer's event was ``pull_request``. Producer runs started
by ``push`` or ``workflow_dispatch`` create a skipped publisher job and never
request a pull request number.

Token responsibilities are deliberately narrow:

.. list-table:: Publisher token responsibilities
   :header-rows: 1
   :widths: 22 28 50

   * - Token
     - Required access
     - Purpose
   * - ``github.token`` (the per-job ``GITHUB_TOKEN``)
     - ``actions: read`` and ``pull-requests: read``
     - PR lookup and artifact download. Some publishers also grant
       ``contents: read`` or ``issues: read`` for metadata.
   * - ``BOT_TOKEN``
     - Write access required by the individual destination
     - Publish previews or benchmark data and create or update the bot comment.
       This token is passed only to the publishing or commenting step and never
       to the producer workflow.

The following pairs use this setup:

.. list-table:: Pull-request producer and publisher workflows
   :header-rows: 1
   :widths: 25 31 44

   * - Producer
     - Publisher
     - Published result
   * - ``docs``
     - ``publish-docs``
     - Documentation preview and status comment
   * - ``benchmarks``
     - ``publish-benchmark-comparison``
     - ASV comparison page and result comment
   * - ``compare-regdata``
     - ``publish-regdata-comparison``
     - Regression-data comparison page and comment
   * - ``codestyle``
     - ``publish-codestyle``
     - Ruff output comment
   * - ``mailmap``
     - ``publish-mailmap``
     - Failure guidance comment
   * - ``orcid-check``
     - ``publish-orcid``
     - Missing-ORCID guidance comment

GitHub displays the producer and publisher as separate workflow runs. Each
publisher's run name includes the pull request title from
``workflow_run.display_title``, and successful bot comments link the source and
publisher runs. GitHub does not provide a token setting that turns a
``workflow_run`` run into a child job of the producer or natively attaches it to
the pull request.

Testing a producer/publisher change
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The producer can be tested in the pull request that changes it. The publisher
is different: GitHub runs the version of a ``workflow_run`` workflow on the
default branch. Consequently, an end-to-end publisher test is only meaningful
after the publisher change is on ``master`` (or in a temporary test repository
whose default branch contains the change).

After the change reaches the default branch:

1. Open or update a pull request from a fork. A fork is important because it
   exercises the empty ``workflow_run.pull_requests`` case.
2. Confirm that the producer workflow checks out and tests the pull-request
   commit successfully.
3. Open the corresponding publisher run. Its title should contain the pull
   request title. The ``Get PR number`` step should identify one current pull
   request and expose a matching head SHA.
4. Confirm that the bot comment links both workflow runs and that any preview
   URL uses the correct pull-request number.
5. Push another commit while an older producer is finishing. The old publisher
   should fail its head-SHA guard and publish nothing; the publisher for the
   newest SHA should publish normally.

When troubleshooting:

- ``403 Resource not accessible by integration`` while looking up a pull request
  usually means ``pull-requests: read`` is missing. The same error while
  downloading an artifact usually means ``actions: read`` is missing.
- ``Artifact not found`` means the producer skipped the upload, used a different
  artifact name, or did not finish successfully. Compare the name in the
  producer with the publisher's name, including the head SHA.

``pull_request_target`` workflows require a separate review. The ``utility``
workflow performs API-only labeling and welcome comments and does not check out
pull-request code. The ``release`` and ``clean-docs`` workflows check out the
trusted default branch; neither selects ``github.event.pull_request.head.sha``.
Never add a pull-request-head checkout to one of these trusted workflows.


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

The documentation preview pipeline builds documentation for pushes to ``master``
and for pull requests. Pull-request previews are published when the
``build-docs`` label is present or when documentation files have changed. Pushes
to ``master`` and manual dispatches also build and publish documentation.

Testing Pipeline
^^^^^^^^^^^^^^^^

The testing pipeline runs one test job per configured operating system:

- Ubuntu tests with the ``not`` continuum marker.
- macOS tests with the ``not`` continuum marker.

When the ``pip-git-tests`` label is present on a pull request, an additional
Git-installed test variant is added for each operating system. The pipeline
includes environment installation, regression data configuration, and coverage
report upload after tests finish.

The ``full-tests`` label enables the separate ``full-tests`` workflow, which runs
the complete test suite and JIT-disabled tests on a self-hosted runner.

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
merged, when invoked by the pre-release workflow, or when manually dispatched.

1. Check out TARDIS with fetch depth 0.
2. Set up Python.
3. Install ``setuptools_scm``.
4. Get the current TARDIS version using ``setuptools_scm`` via a helper script.
5. Get the next TARDIS version using ``setuptools_scm``.
6. Create a GitHub release using the new version as the tag.
7. Wait 2 minutes for Zenodo to update the new TARDIS release.
8. Fetch the new DOI from Zenodo using the Zenodo API and create a badge.
9. Generate the changelog with ``orhun/git-cliff-action`` using the TARDIS
   ``pyproject.toml`` configuration.
10. Update the release description with the changelog and Zenodo badge.
11. Include environment lock files in release assets.

Post-Release
^^^^^^^^^^^^

The post-release action updates the changelog, citation, and credits in the main
repository.

Changelog job:

1. Check out TARDIS with fetch depth 0.
2. Get the current release tag.
3. Generate a changelog with ``orhun/git-cliff-action``.
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
7. Dispatch updates to the TARDIS website with
   ``peter-evans/repository-dispatch``.

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

Regression Data Comparison Workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Regression Data Comparison workflow compares regression data between the
current branch and the base branch on pull requests. It runs only on pull
requests, not on ``master``.

The workflow generates regression data for the latest commit on the pull request
and compares it with ``master`` using the comparison notebook. The notebook is
uploaded as an artifact and pushed to the ``reg-data-comp`` repository for
previews in the bot comment.

The comparison uses the :ref:`secure-pull-request-publishing` split.
``compare-regdata`` runs pull-request code without repository secrets and
uploads the comparison artifacts. ``publish-regdata-comparison`` then checks
the pull-request label and publishes the artifacts and bot comment.

The workflow exports images from the comparison notebook and embeds them in the
bot comment. Unless there are key changes to HDF files in the regression data,
the bot shows two images: spectrum change and relative changes in keys. If
there are key changes, the bot shows a third image visualizing key changes.

LFS-Cache Workflow
^^^^^^^^^^^^^^^^^^

The ``LFS-cache`` workflow caches regression data and atomic data. It can be
triggered manually or called by another workflow. It checks for a matching
cache, performs an LFS pull when the cache is missing, and saves the resulting
objects. The ``setup_lfs`` action restores cached objects for consuming jobs and
fails when the required cache is unavailable.

Research Papers Workflow
^^^^^^^^^^^^^^^^^^^^^^^^

The ``tardis-research-papers`` workflow runs monthly and can also be dispatched
manually. It:

1. Runs the ADS notebook with the repository's research-paper environment.
2. Creates a pull request with updated research-paper data.
3. Enables auto-merge after the required checks pass.
4. Dispatches a ``fetch-papers`` event to ``tardis-sn/tardis-org-data``.
