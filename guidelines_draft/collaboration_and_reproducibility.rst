.. _collaboration_and_reproducibility:

Collaboration And Reproducibility
=================================

Guidance for issues, Git workflow, pull requests, CI, releases, and collaboration infrastructure.

.. _continuous-integration:

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

.. raw:: html

   <span style="color:red">Deleted: removed repeated GitHub Actions overview sentence already covered in the explanation.</span>


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


.. raw:: html

   <span style="color:red">Added: example from a real TARDIS workflow.</span>


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

.. _developer-workflow:

Developer Workflow
------------------

.. _tutorial-add-a-new-feature:

Tutorial: Add a New Feature
~~~~~~~~~~~~~~~~~~~~~~~~~~~

TARDIS development aims to follow test-driven development. A new feature should
move from a documented need, to a failing test, to implementation, to passing
tests, to documentation and changelog updates.

1. Document the feature to be added in an issue.
2. If you are unsure whether the feature already exists, ask the mailing list.
3. Write a test that demonstrates the feature that will be added.
4. Run the test to verify that it fails in the way you expect.
5. If the test fails unexpectedly, the test may be wrong. Ask the group for
   guidance.
6. If the test passes, the feature may already exist. In that case, you have
   added test coverage for existing behavior.
7. Add the feature.
8. Run the test to verify that it passes.
9. Write documentation about the feature.
10. Close the issue or partial pull request and add the change to the changelog.

Large features should usually be broken into small, quantifiable goals that can
be acted on separately.

.. raw:: html

   <span style="color:red">Added: concrete feature-development example.</span>


For example, if a change adds support for a new density configuration, a small
development path would be: open an issue describing the new density option, add
a failing test near the existing model or configuration tests, update
``tardis/io/model/parse_density_configuration.py``, run the targeted tests, and
then document the new configuration behavior.

.. _how-to-guide-introduce-a-new-feature:

How-To Guide: Introduce a New Feature
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before implementing a new feature, make the feature reviewable:

1. Open or find an issue that states the user need, the expected behavior, and
   any scientific or data assumptions.
2. Identify the smallest code path that owns the behavior. For example, a new
   model input format likely belongs near ``tardis/io/model/``, while a new
   visualization behavior likely belongs near ``tardis/visualization/``.
3. Add a failing unit test or regression test that captures the new behavior.
4. Implement the smallest code change that makes the test pass.
5. Add or update documentation for the feature.
6. Run targeted tests, Ruff, and any relevant documentation build.
7. Open a pull request with the issue link, test commands, and reviewer notes.

Example feature checklist:

.. code-block:: markdown

   Issue
   - Adds support for a new density option in CSVY model input.

   Code
   - Update tardis/io/model/parse_density_configuration.py.
   - Add or update tests under tardis/io/model/readers/tests.

   Validation
   - pytest tardis/io/model/readers/tests
   - ruff check tardis/io/model

   Documentation
   - Update the relevant configuration or model reader page.


.. raw:: html

   <span style="color:red">Added: explicit feature-introduction workflow.</span>


.. _explanation-developer-workflow:

Explanation: Developer Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The developer workflow is intended for anyone contributing to TARDIS
development. Much of the workflow is taken from Astropy's development workflow,
with credit to the Astropy team.

The workflow emphasizes:

- Filing issues for bugs, inconsistencies, missing functionality, and planned
  features.
- Working from a fork.
- Keeping an ``upstream`` remote connected to the main TARDIS repository.
- Avoiding development on local ``master``.
- Starting each separable task from current upstream trunk.
- Keeping work organized through feature branches.
- Running tests before review.
- Updating documentation.
- Opening pull requests for review and merge.

This helps maintain readable history and helps maintainers understand what
changed and why.

.. _explanation-test-driven-development:

Explanation: Test-Driven Development
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TARDIS aims to use test-driven development. The preferred path starts by
describing a feature, then writing a test for the missing behavior, then
implementing the feature until the test passes.

Unexpected test failures can indicate that the test is wrong. A passing test
before implementation can mean the feature already exists and only lacked
coverage.

.. _external-links:

External Links
--------------

.. _reference-external-link-reference:

Reference: External Link Reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Project links:

- TARDIS: http://tardis.readthedocs.org
- TARDIS GitHub: https://github.com/tardis-sn/tardis
- TARDIS developer mailing list:
  https://groups.google.com/forum/#!forum/tardis-sn-dev

Git and workflow links:

- Git: http://git-scm.com/
- GitHub: http://github.com
- GitHub Help: http://help.github.com
- GitHub fork help: https://help.github.com/fork-a-repo/
- Git SSH keys: https://help.github.com/articles/generating-ssh-keys
- Git cheat sheet: http://cheat.errtheblog.com/s/git
- Pro Git book: http://progit.org/
- Git tutorial: http://schacon.github.com/git/gittutorial.html
- Git user manual: http://schacon.github.com/git/user-manual.html
- Rebase without tears:
  http://matthew-brett.github.com/pydagogue/rebase_without_tears.html
- Resolving a merge:
  http://schacon.github.com/git/user-manual.html#resolving-a-merge

Tool links:

- Python: http://www.python.org
- Conda: https://docs.conda.io/en/latest/
- Miniconda: https://docs.anaconda.com/free/miniconda/
- Mamba: https://github.com/mamba-org/mamba
- Mini-forge: https://github.com/conda-forge/miniforge
- Pytest: https://docs.pytest.org/en/latest/
- Ruff: https://docs.astral.sh/ruff/
- Pre-commit: https://pre-commit.com/
- Sphinx: https://www.sphinx-doc.org/
- Sphinx quickstart:
  https://www.sphinx-doc.org/en/master/usage/quickstart.html
- Sphinx reStructuredText basics:
  https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html
- Sphinx toctree directive:
  https://www.sphinx-doc.org/en/master/usage/restructuredtext/directives.html#directive-toctree
- Sphinx build options:
  https://www.sphinx-doc.org/en/master/man/sphinx-build.html#options
- nbsphinx: https://nbsphinx.readthedocs.io/
- nbsphinx hidden cells:
  https://nbsphinx.readthedocs.io/en/0.8.7/hidden-cells.html
- Jupyter: https://jupyter.org/
- Binder: https://mybinder.org/
- NumPy docstring format:
  https://numpydoc.readthedocs.io/en/latest/format.html
- AirSpeed Velocity: https://asv.readthedocs.io/en/latest/index.html
- ASV installation:
  https://asv.readthedocs.io/en/latest/installing.html#installing-airspeed-velocity

.. raw:: html

   <span style="color:red">Deleted: removed several older general-purpose Git learning links from this draft to keep the reference focused on links directly used by the developer workflow.</span>


TARDIS data and CI links:

- TARDIS regression data:
  https://github.com/tardis-sn/tardis-regression-data
- Git LFS tutorial:
  https://www.atlassian.com/git/tutorials/git-lfs
- Tests workflow:
  https://github.com/tardis-sn/tardis/blob/master/.github/workflows/tests.yml
- GitHub Actions:
  https://docs.github.com/en/actions
- GitHub-hosted runners:
  https://docs.github.com/en/actions/using-github-hosted-runners/about-github-hosted-runners/about-github-hosted-runners#standard-github-hosted-runners-for-public-repositories
- Pull request event:
  https://docs.github.com/en/actions/using-workflows/events-that-trigger-workflows#pull_request
- Push event:
  https://docs.github.com/en/actions/using-workflows/events-that-trigger-workflows#push
- Carsus: https://tardis-sn.github.io/carsus/

Matterbridge links:

- OpenSupernova.org: http://opensupernova.org
- Private Matterbridge repository:
  https://github.com/tardis-sn/tardis-matterbridge
- Upstart script:
  https://www.digitalocean.com/community/tutorials/the-upstart-event-system-what-it-is-and-how-to-use-it
- Systemd service:
  https://freedesktop.org/software/systemd/man/systemd.service.html
- Slack setup:
  https://github.com/42wim/matterbridge/wiki/Slack-bot-setup
- Mattermost setup:
  https://github.com/42wim/matterbridge/wiki/Section-Mattermost-%28basic%29
- Slack scopes comment:
  https://github.com/42wim/matterbridge/issues/964#issuecomment-612721850
- Gitter authentication:
  https://developer.gitter.im/docs/authentication
- Matterbridge releases:
  https://github.com/42wim/matterbridge/releases/latest

.. _git-workflow:

Git Workflow
------------

.. _how-to-guide-connect-your-fork-to-upstream:

How-To Guide: Connect Your Fork to Upstream
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Add the main TARDIS repository as ``upstream``:

.. code-block:: shell

   git remote add upstream https://github.com/tardis-sn/tardis.git


The ``upstream`` name is the convention used for the main TARDIS repository. The
workflow uses the HTTPS URL for upstream because it is read-only by default for
most contributors and avoids accidental writes to the main repository.

Confirm the remotes:

.. code-block:: shell

   git remote -v show


Expected output should include:

.. code-block:: text

   upstream   https://github.com/tardis-sn/tardis.git (fetch)
   upstream   https://github.com/tardis-sn/tardis.git (push)
   origin     git@github.com:your-user-name/tardis.git (fetch)
   origin     git@github.com:your-user-name/tardis.git (push)


.. _how-to-guide-fork-and-clone-tardis:

How-To Guide: Fork and Clone TARDIS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You only need to fork once for each package you want to contribute to.

1. Log into GitHub.
2. Go to the TARDIS GitHub home page.
3. Click the fork button.
4. Clone your fork:

   .. code-block:: shell

      git clone git@github.com:your-user-name/tardis.git
      cd tardis


5. Check branches:

   .. code-block:: shell

      git branch -a


6. Check remotes:

   .. code-block:: shell

      git remote -v


At this point, ``origin`` should point to your GitHub fork.

.. _explanation-git-branching-strategy:

Explanation: Git Branching Strategy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the developer workflow, ``master`` is referred to as the trunk. Feature work
should begin from ``upstream/master``.

Feature branches make reviews easier because each branch contains a related set
of edits. Contributors should avoid merging trunk or other branches into a
feature branch when possible. If trunk changes need to be incorporated, rebasing
is preferred because it replays feature commits on top of the latest trunk and
keeps history easier to read.

A rebase transforms a history like:

.. code-block:: text

         A---B---C cool-feature
        /
   D---E---F---G trunk


into:

.. code-block:: text

                 A'--B'--C' cool-feature
                /
   D---E---F---G trunk


Merge conflicts can occur if the feature branch and trunk changed the same
files. Resolve conflicts using the Git rebase documentation and related merge
resolution guidance.

.. _reference-git-command-reference:

Reference: Git Command Reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Clone a fork:

.. code-block:: shell

   git clone git@github.com:your-user-name/tardis.git


Inspect branches:

.. code-block:: shell

   git branch -a


Inspect remotes:

.. code-block:: shell

   git remote -v
   git remote -v show


Add upstream:

.. code-block:: shell

   git remote add upstream https://github.com/tardis-sn/tardis.git


Fetch upstream:

.. code-block:: shell

   git fetch upstream


Start a branch from trunk:

.. code-block:: shell

   git checkout upstream/master
   git checkout -b my-new-feature


Push a branch:

.. code-block:: shell

   git push origin my-new-feature
   git push --set-upstream origin my-new-feature


Check status and diff:

.. code-block:: shell

   git status
   git diff


Stage files:

.. code-block:: shell

   git add new_file_name
   git add modified_file_name


Commit:

.. code-block:: shell

   git commit -m "A commit message"


Rebase:

.. code-block:: shell

   git fetch upstream
   git checkout cool-feature
   git branch tmp cool-feature
   git rebase upstream/master


Abort a rebase:

.. code-block:: shell

   git rebase --abort


Force push a rebased branch to your fork:

.. code-block:: shell

   git push -f origin cool-feature


Reset to a backup branch:

.. code-block:: shell

   git reset --hard tmp


Inspect reflog:

.. code-block:: shell

   git reflog show cool-feature


.. _tutorial-prepare-your-local-development-fork:

Tutorial: Prepare Your Local Development Fork
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This path prepares a new contributor to work on TARDIS locally.

1. Set up a Python environment. The original developer workflow recommends
   Anaconda and refers readers to the installation guide.
2. Create a GitHub account if you do not already have one.
3. Configure your GitHub account for write access, including SSH keys.
4. Fork the TARDIS repository from the TARDIS GitHub home page.
5. Clone your fork:

   .. code-block:: shell

      git clone git@github.com:your-user-name/tardis.git
      cd tardis


6. Inspect local and remote branches:

   .. code-block:: shell

      git branch -a
      git remote -v


7. Add the main TARDIS repository as ``upstream``:

   .. code-block:: shell

      git remote add upstream https://github.com/tardis-sn/tardis.git


8. Confirm the remote configuration:

   .. code-block:: shell

      git remote -v show


   Expected remotes include ``origin``, pointing to your fork, and ``upstream``,
   pointing to ``https://github.com/tardis-sn/tardis.git``.

9. Install TARDIS in development mode:

   .. code-block:: shell

      pip install -e .


   This installs TARDIS so imports use your repository clone regardless of your
   working directory. Edits in your clone are available the next time you start a
   Python interpreter and ``import tardis``.

.. _how-to-guide-recover-from-git-mistakes:

How-To Guide: Recover from Git Mistakes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a rebase goes wrong while it is in progress:

.. code-block:: shell

   git rebase --abort


If you notice the problem after the rebase and made a backup branch:

.. code-block:: shell

   git reset --hard tmp


If you forgot to make a backup branch, inspect the reflog:

.. code-block:: shell

   git reflog show cool-feature


Example reflog:

.. code-block:: text

   8630830 cool-feature@{0}: commit: BUG: io: close file handles immediately
   278dd2a cool-feature@{1}: rebase finished: refs/heads/my-feature-branch onto 11ee694744f2552d
   26aa21a cool-feature@{2}: commit: BUG: lib: make seek_gzip_factory not leak gzip obj


Reset to the point before the bad rebase:

.. code-block:: shell

   git reset --hard cool-feature@{2}


.. _how-to-guide-set-up-github-access:

How-To Guide: Set Up GitHub Access
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Create a GitHub account at https://github.com if needed.
2. Configure write access with SSH keys.
3. Use GitHub Help's SSH key instructions:
   https://help.github.com/articles/generating-ssh-keys

.. _how-to-guide-start-a-feature-branch:

How-To Guide: Start a Feature Branch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the developer workflow, the TARDIS ``master`` branch is the trunk.

1. Do not use your local ``master`` branch for development. Consider deleting it
   to reduce confusion.
2. Fetch the latest upstream changes:

   .. code-block:: shell

      git fetch upstream


3. Start from the current upstream trunk:

   .. code-block:: shell

      git checkout upstream/master


4. Create a new feature branch:

   .. code-block:: shell

      git checkout -b my-new-feature


Use a new branch for each separable set of changes: one task, one branch. Choose
an informative name, such as ``bugfix-for-issue-14``,
``refactor-density-parser``, or ``update-regression-data-docs``.

.. raw:: html

   <span style="color:red">Added: replaced generic branch-name examples with TARDIS-oriented examples.</span>


Push the branch to your fork:

.. code-block:: shell

   git push origin my-new-feature


With Git 1.7 or newer, you can set the upstream branch:

.. code-block:: shell

   git push --set-upstream origin my-new-feature


.. _how-to-guide-use-the-editing-workflow:

How-To Guide: Use the Editing Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Make your changes.
2. Run tests to check for regressions:

   .. code-block:: shell

      pytest tardis


3. If Sphinx is installed, check the documentation build:

   .. code-block:: shell

      cd docs
      make html


   The build should succeed and should not report warnings.

4. Check changed files:

   .. code-block:: shell

      git status


5. Inspect the actual changes:

   .. code-block:: shell

      git diff


6. Add new files:

   .. code-block:: shell

      git add new_file_name


7. Add modified files you want to commit:

   .. code-block:: shell

      git add modified_file_name


8. Check what will be committed:

   .. code-block:: shell

      git status


9. Commit:

   .. code-block:: shell

      git commit -m "A commit message"


10. Push to your fork:

    .. code-block:: shell

       git push


.. raw:: html

   <span style="color:red">Added: concrete editing workflow example.</span>


For a small parser change, a typical local loop might be:

.. code-block:: shell

   ruff check tardis/io/model/parse_density_configuration.py
   pytest tardis/io/model/readers/tests
   git diff
   git add tardis/io/model/parse_density_configuration.py
   git commit -m "Improve density parser validation"
   git push

.. _issues:


Issues
------

.. _how-to-guide-report-issues:

How-To Guide: Report Issues
~~~~~~~~~~~~~~~~~~~~~~~~~~~

TARDIS is under constant development, so bugs, inconsistencies, and missing
functionality can occur. File issues on the official GitHub repository:

https://github.com/tardis-sn/tardis

Interested users are encouraged to contribute to TARDIS development by following
the development workflow.

.. _matterbridge:


Matterbridge
------------

.. _how-to-guide-edit-the-cron-job:

How-To Guide: Edit the Cron Job
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SSH to OpenSupernova.org and run:

.. code-block:: shell

   sudo crontab -e


.. _explanation-matterbridge:

Explanation: Matterbridge
~~~~~~~~~~~~~~~~~~~~~~~~~

Matterbridge connects messaging channels across platforms such as Slack,
Mattermost, and Gitter. The ``matterbridge`` binary is used with a TOML
configuration file:

.. code-block:: shell

   ./matterbridge -conf config-tardis-matterbridge.toml


The TOML file contains parameters required to connect rooms, including tokens
and passwords. When the application runs, messages can be shared between the
connected rooms.

TARDIS keeps a service running on the OpenSupernova.org server to run
Matterbridge as a daemon. Configuration files are stored in a private GitHub
repository, including the custom TOML and Linux service files.

The server runs Ubuntu 14.04, so TARDIS uses an Upstart script instead of a
Systemd service file. A Systemd file is also included in the repository. A cron
job restarts the service periodically to prevent disconnections.

.. _reference-matterbridge-reference:

Reference: Matterbridge Reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Matterbridge command:

.. code-block:: shell

   ./matterbridge -conf config-tardis-matterbridge.toml


Matterbridge configuration:

- Use the TOML file in ``tardis-matterbridge`` as an example.
- Configure as many gateways as needed.
- For Slack, follow the Matterbridge wiki Slack setup steps and read the linked
  comment warning not to add scopes manually.
- For Gitter, create a dedicated GitHub bot account and copy the token.
- For Mattermost, follow the Matterbridge wiki Mattermost setup steps.

Important paths:

- Matterbridge executable: ``/usr/local/bin/matterbridge``
- Matterbridge TOML config:
  ``/usr/local/etc/matterbridge/config-tardis-matterbridge.toml``
- Upstart service config: ``/etc/init/matterbridge.conf``

Important repositories and services:

- OpenSupernova.org server: http://opensupernova.org
- Private Matterbridge repository:
  https://github.com/tardis-sn/tardis-matterbridge
- Matterbridge releases: https://github.com/42wim/matterbridge/releases/latest

.. _how-to-guide-set-up-matterbridge:

How-To Guide: Set Up Matterbridge
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. raw:: html

   <span style="color:red">Deleted: removed repeated Matterbridge overview because the explanation above already defines the service.</span>


1. SSH to the OpenSupernova.org server.
2. Download the Matterbridge binary for Linux from the releases section.
3. Make the file executable and rename it to ``matterbridge``.
4. Copy the ``matterbridge`` executable to ``/usr/local/bin``.
5. Clone the ``tardis-matterbridge`` repository in ``$HOME``.
6. Copy ``config-tardis-matterbridge.toml`` to
   ``/usr/local/etc/matterbridge/``.
7. Copy ``matterbridge.conf`` to ``/etc/init/``.
8. Start the service:

   .. code-block:: shell

      sudo service matterbridge start


9. Test the gateways.

.. _how-to-guide-update-server-configuration:

How-To Guide: Update Server Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After updating the TOML file:

1. SSH to the OpenSupernova.org server.
2. Copy the new ``config-tardis-matterbridge.toml`` to
   ``/usr/local/etc/matterbridge/``.
3. Restart the service:

   .. code-block:: shell

      sudo service matterbridge restart


4. Test the gateways.
5. If everything works, open a pull request to ``tardis-matterbridge`` with the
   new TOML file.

.. _pull-requests:

Pull Requests
-------------

.. _how-to-guide-open-a-pull-request:

How-To Guide: Open a Pull Request
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When you are ready for review or merge consideration:

1. Go to your fork, for example ``https://github.com/your-user-name/tardis``.
2. Use the branch dropdown to select the branch with your changes.
3. Click the pull request button.
4. Enter a title and an explanation of what you have done.
5. Include anything that needs particular attention, such as a complicated
   change or code you are unsure about.
6. If the request is not ready to merge, say so in the pull request message.
   Opening an unfinished pull request can still be a good way to start a
   preliminary code review.
7. Build documentation for the pull request using the documentation preview
   process so reviewers can see how notebooks, docstrings, and API docs render.

.. raw:: html

   <span style="color:red">Added: example pull request description template.</span>


Example pull request description:

.. code-block:: markdown

   Summary
   - Update density parsing validation.
   - Add a regression or unit test for the new behavior.

   Checks
   - pytest tardis/io/model/readers/tests
   - ruff check tardis/io/model

   Notes for reviewers
   - Please check whether the error message is clear for invalid density types.


.. _how-to-guide-review-pull-requests:

How-To Guide: Review Pull Requests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GitHub supports line comments and discussion threads on pull requests. Sometimes
it is easier to demonstrate a suggested change by making a pull request on a
pull request.

Add pull request refs to the ``upstream`` remote in ``.git/config``:

.. code-block:: ini

   [remote "upstream"]
       url = git@github.com:tardis-sn/tardis.git
       fetch = +refs/heads/*:refs/remotes/upstream/*
       fetch = +refs/pull/*/head:refs/remotes/upstream/pr/*


Fetch upstream:

.. code-block:: shell

   git fetch upstream


You can then check out a pull request branch, for example pull request 116:

.. code-block:: shell

   git checkout upstream/pr/116


This leaves you in a detached ``HEAD`` state. Create a branch for your changes:

.. code-block:: shell

   git checkout -b helping-with-PR116


After making and committing your changes, push them to your fork:

.. code-block:: shell

   git push origin helping-with-PR116


.. _how-to-guide-update-a-pull-request:

How-To Guide: Update a Pull Request
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pull requests often remain open while feedback is added and other pull requests
are merged. If changes on trunk affect your work, rebase your feature branch on
top of the current trunk.

Fetch trunk and switch to your branch:

.. code-block:: shell

   git fetch upstream
   git checkout cool-feature


Create a temporary backup branch:

.. code-block:: shell

   git branch tmp cool-feature


Rebase onto upstream trunk:

.. code-block:: shell

   git rebase --onto upstream/master upstream/master cool-feature


If you are already on ``cool-feature``, the command can be shortened:

.. code-block:: shell

   git rebase upstream/master


When everything looks good, delete the backup branch:

.. code-block:: shell

   git branch -D tmp


If the feature branch is already on GitHub, force push after rebasing:

.. code-block:: shell

   git push -f origin cool-feature


Force pushing overwrites the branch on GitHub and can lose commits if used
incorrectly. Never force push to the main TARDIS repository, typically called
``upstream``, because that rewrites shared history.
