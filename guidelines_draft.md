# TARDIS Developer Documentation Draft

## Table of Contents

- [Code Style And Syntax](#code-style-and-syntax)
  - [Code Quality](#code-quality)
  - [Naming Conventions](#naming-conventions)
  - [Pre-Commit](#pre-commit)
  - [Ruff](#ruff)
  - [Typed NumPy Arrays](#typed-numpy-arrays)
- [Collaboration And Reproducibility](#collaboration-and-reproducibility)
  - [Continuous Integration](#continuous-integration)
  - [Developer Workflow](#developer-workflow)
  - [External Links](#external-links)
  - [Git Workflow](#git-workflow)
  - [Issues](#issues)
  - [Matterbridge](#matterbridge)
  - [Pull Requests](#pull-requests)
- [Data And Codebase Structure](#data-and-codebase-structure)
  - [Developer FAQ](#developer-faq)
  - [Development Install](#development-install)
  - [Numba Debugging](#numba-debugging)
- [Documentation](#documentation)
  - [Docstrings](#docstrings)
  - [Documentation Builds](#documentation-builds)
  - [Documentation Pages](#documentation-pages)
- [Testing And Validation](#testing-and-validation)
  - [Benchmarks](#benchmarks)
  - [Regression Data](#regression-data)
  - [Test Coverage](#test-coverage)
  - [Tests](#tests)

## Code Style And Syntax

Guidance for writing readable, maintainable TARDIS code.

### Code Quality

#### Explanation: Code Quality

Code quality helps new developers understand previous work. In open-source
software such as TARDIS, high-quality code is essential.

High-quality code:

- Does what it is supposed to do.
- Does not contain defects or problems.
- Handles edge cases safely.
- Raises or handles exceptions deliberately.
- Is easy to read, maintain, and extend.

TARDIS follows PEP 8. Ruff handles many style concerns such as whitespace,
string quotes, and code layout, but developers are still responsible for naming
conventions:

- Function names should be lowercase with words separated by underscores.
- Variable names use the same convention as function names.
- Class names should use CapWords.

<span style="color:red">Added: concrete TARDIS naming example.</span>

For example, `tardis/io/model/parse_density_configuration.py` uses descriptive
function names such as `parse_density_section_config`,
`parse_density_from_csvy`, and `calculate_power_law_density`. These names are
preferable to short generic names because they state the model component and
operation being performed.

##### Edge Cases and Exceptions

Code should anticipate errors that may occur during execution. If an exception
is likely and can be handled, handle it. If an edge case would cause incorrect
behavior, raise an appropriate exception with a useful message.

Example:

```python
def _calculate_plotting_data(self, packets_mode, packet_wvl_range, distance):
    if packets_mode not in ["virtual", "real"]:
        raise ValueError(
            "Invalid value passed to packets_mode. Only "
            "allowed values are 'virtual' or 'real'"
        )
    # Rest of the code ...
```

Here, `packets_mode` must be either `"virtual"` or `"real"`. A specific
`ValueError` tells the user what went wrong and prevents invalid code from
continuing.

#### Reference: Code Quality Reference

TARDIS follows PEP 8.

Rule locations:

- Permanent Ruff rules: `pyproject.toml`
- Non-permanent Ruff rules: `.ruff.toml`

Naming conventions:

- Functions: lowercase with underscores.
- Variables: lowercase with underscores.
- Classes: CapWords.

Pre-commit install:

```shell
pip install pre-commit
pre-commit install
```

<span style="color:red">Deleted: removed duplicate Ruff command examples from
the reference because the Ruff how-to gives runnable commands.</span>

### Naming Conventions

#### Explanation: Variable and Function Names

Use names that describe the scientific or codebase concept, not only the local
implementation detail. Prefer names such as `density_configuration`,
`time_explosion`, `velocity_field_index`, and `density_field_index` when those
objects refer to specific model quantities or schema locations.

Avoid names that are too generic for review, such as `data`, `arr`, `x`, or
`idx`, unless the scope is very small and the meaning is obvious from the
surrounding line. If a variable stores the position of a named field, use the
`_index` suffix and include the concept being indexed:

```python
velocity_field_index = [
    field.name for field in csvy_model_config.datatype.fields
].index("velocity")

density_field_index = [
    field.name for field in csvy_model_config.datatype.fields
].index("density")
```

Use `_idx` only for short-lived loop or array positions where the indexed object
is already named nearby. Use `_index` for persistent variables, pandas indexes,
schema indexes, or public-facing names where readability matters more than
brevity.

Good examples:

```python
velocity_field_index = fields.index("velocity")
density_field_index = fields.index("density")
packet_idx = packet_indices[i]
```

Poor examples:

```python
idx = fields.index("density")
i = packet_indices[i]
data = parse_density_section_config(...)
```

<span style="color:red">Added: naming convention guidance for `_idx`,
`_index`, and TARDIS-style variable/function names.</span>

### Pre-Commit

#### How-To Guide: Use Pre-Commit

Pre-commit hooks are optional local tools that run checks before commits.

Install pre-commit:

```shell
pip install pre-commit
```

Set up hooks for the repository:

```shell
pre-commit install
```

This only needs to be done once per repository. The hooks then run automatically
on each commit.

<span style="color:red">Added: example command for checking a TARDIS file before
committing.</span>

To run the hooks manually on a specific file, use:

```shell
pre-commit run --files tardis/io/model/parse_density_configuration.py
```

### Ruff

#### How-To Guide: Run Ruff

TARDIS follows PEP 8 and uses Ruff for linting and formatting.

Install Ruff:

```shell
conda install -c conda-forge ruff
```

Lint code:

```shell
ruff check <source_file_or_directory>
```

Lint and fix automatically fixable issues:

```shell
ruff check <source_file_or_directory> --fix
```

<span style="color:red">Added: TARDIS-specific Ruff examples.</span>

For example, lint one parser module:

```shell
ruff check tardis/io/model/parse_density_configuration.py
```

Or lint and fix a package area:

```shell
ruff check tardis/io/model --fix
```

TARDIS adopts linting rules used by Astropy. Permanent rules are defined in
`pyproject.toml`. Non-permanent rules are defined in `.ruff.toml`. Add new rules
to `.ruff.toml`. To add a permanent rule, open a pull request against
`pyproject.toml`.

### Typed NumPy Arrays

#### Reference: NumPy Array Typing

Use `numpy.typing.NDArray` for typed NumPy array annotations in new or touched
code when the function accepts or returns concrete NumPy arrays. NumPy provides
`numpy.typing.NDArray` as a runtime-available alias for annotating arrays with a
dtype and unspecified shape.

```python
from typing import Any

import numpy as np
import numpy.typing as npt


def normalize(values: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    return values / np.sum(values)


def as_array(values: npt.ArrayLike) -> npt.NDArray[Any]:
    return np.asarray(values)
```

Prefer typed arrays over bare `np.ndarray` when the dtype is meaningful to the
reader or reviewer:

```python
def calculate_density_after_time(
    density_0: npt.NDArray[np.float64],
    time_0: float,
    time_explosion: float,
) -> npt.NDArray[np.float64]:
    ...
```

The typed NumPy API is stricter than runtime NumPy. For example, type checkers
will discourage patterns that create object arrays accidentally or mutate array
dtypes directly. If a function accepts broad array-like input, use
`npt.ArrayLike` for the input and return a typed `npt.NDArray[...]` after
conversion.

See the official NumPy typing reference:
https://numpy.org/doc/stable/reference/typing.html

<span style="color:red">Added: typed NumPy ndarray guidance and examples.</span>

## Collaboration And Reproducibility

Guidance for issues, Git workflow, pull requests, CI, releases, and collaboration infrastructure.

### Continuous Integration

#### Explanation: Continuous Integration

TARDIS uses continuous integration. When a change is proposed by pull request or
merged into `master`, a service clones the repository, checks out the current
commit, and runs TARDIS tests. This helps detect bugs immediately.

TARDIS currently uses GitHub Actions for pipelines. A workflow is a YAML
configuration file with sections such as variables, jobs, and steps. Workflows
run commands when triggered by events, such as pushes or pull requests.

Making changes to an existing pipeline is done through a pull request. Creating
a new workflow requires adding a YAML file under `.github/workflows`.

#### How-To Guide: Configure a GitHub Actions Workflow

To create a workflow, add a YAML file under `.github/workflows`.

<span style="color:red">Deleted: removed repeated GitHub Actions overview
sentence already covered in the explanation.</span>

Use common TARDIS setup actions and settings:

1. Use the `setup_lfs` action and `lfs-cache` workflow when regression or atomic
   data is needed.
2. Use the `lfs-cache` workflow to cache regression data and atomic data and to
   check whether the cache is available.
3. Use the `setup_env` action to configure variables and settings for the
   pipeline.
4. Set the shell in the YAML configuration:

   ```yaml
   defaults:
     run:
       shell: bash -l {0}
   ```

<span style="color:red">Added: example from a real TARDIS workflow.</span>

The `docs` workflow in `.github/workflows/build-docs.yml` uses the shared LFS
cache workflow, restores sparse atomic data, sets up the environment, and builds
the docs:

```yaml
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
```

#### Reference: Continuous Integration Reference

##### Cache Keys

Regression data cache key format:

```text
tardis-regression-<data-type>-<hash>-v1
```

Examples:

- `tardis-regression-atom-data-sparse-<hash>-v1`: atomic data cache.
- `tardis-regression-full-data-<hash>-v1`: full TARDIS regression data cache.

Used in the `setup_lfs` action.

Environment cache key format:

```text
tardis-conda-env-<os-label>-<hash>-v1
```

Examples:

- `tardis-conda-env-linux-<hash>-v1`: Linux conda environment.
- `tardis-conda-env-macos-<hash>-v1`: macOS conda environment.

Used in the `setup_env` action.

The `-v1` suffix allows future cache invalidation. The `lfs-cache` workflow
fails if the cache is unavailable and does not pull LFS data by default. If the
`allow_lfs_pull` label is added to a pull request, the workflow pulls LFS data.
Use this label sparingly and only with caution.

##### Common Pipeline Steps

- Use `setup_lfs` and the `lfs-cache` workflow when regression or atomic data is
  required.
- Use `setup_env` to configure environment variables and settings.
- Use `bash -l {0}` as the run shell in workflow YAML.

##### Documentation Build Pipeline

The documentation build pipeline builds and deploys the TARDIS documentation
website.

##### Documentation Preview Pipeline

The documentation preview pipeline does not run on the main repository, only on
forks. It supports pull request documentation previews.

##### Testing Pipeline

The testing pipeline runs four concurrent jobs:

- Ubuntu tests without the continuum marker.
- Ubuntu tests with the continuum marker.
- macOS tests without the continuum marker.
- macOS tests with the continuum marker.

This is two platforms multiplied by two test types. The pipeline includes
environment installation, regression data configuration, and coverage report
upload after tests finish.

GPU tests can be triggered by applying the `full-tests` label to a pull request.

##### Release Pipeline

TARDIS publishes a new release every Sunday at 00:00 UTC.

##### Pre-Release

The pre-release action clones `tardis-sn/tardis_zenodo`, runs a notebook to
generate `.zenodo.json`, and pushes the file to the root of the TARDIS
repository. The file creates a new TARDIS version on Zenodo with all committers
as authors. A pull request is created and automatically merged if required
checks pass.

Zenodo job:

1. Check out `tardis-sn/tardis_zenodo`.
2. Wait for the Zenodo webhook to be available with a 3 minute sleep.
3. Set up the Python environment stored in `tardis-sn/tardis_zenodo`.
4. Store the Zenodo API secret key in an environment variable.
5. Run the notebook to generate `.zenodo.json`.
6. Re-run if there are errors and ignore errors.
7. Upload `.zenodo.json` as an artifact.

Pip tests job:

Runs the TARDIS test suite after installing TARDIS from `master` with pip using
Git, instead of an editable install. This catches cases where a new module lacks
an `__init__.py` file or data files are missing from
`[tool.setuptools.package-data]` in `pyproject.toml`.

Pre-release pull request job:

1. Relies on Zenodo and pip test steps completing.
2. Checks out the TARDIS repository.
3. Downloads artifacts from previous steps.
4. Checks for `.zenodo.json` and uses it if generated.
5. Gets the current date.
6. Creates a bot pull request on the `tardis-bot` fork using branch
   `pre-release-<date>` with the new `.zenodo.json`.
7. Waits 1 minute for the pull request to be created.
8. Automatically approves the pull request using tokens from infrastructure and
   core coordinator members.
9. Enables auto-merge.

##### Release

The release job creates a GitHub release after the pre-release pull request is
merged.

1. Check out TARDIS with fetch depth 0.
2. Set up Python.
3. Install `setuptools_scm` and `git-cliff`.
4. Get the current TARDIS version using `setuptools_scm` via a helper script.
5. Get the next TARDIS version using `setuptools_scm`.
6. Create a GitHub release using the new version as the tag.
7. Wait 2 minutes for Zenodo to update the new TARDIS release.
8. Fetch the new DOI from Zenodo using the Zenodo API and create a badge.
9. Generate the changelog with `git-cliff`.
10. Update the release description with the changelog and Zenodo badge.
11. Include environment lock files in release assets.

##### Post-Release

The post-release action updates the changelog, citation, and credits in the main
repository.

Changelog job:

1. Check out TARDIS with fetch depth 0.
2. Get the current release tag.
3. Generate a changelog with `git-cliff`.
4. Upload `CHANGELOG.md` as an artifact.

Citation job:

1. Check out TARDIS.
2. Wait 3 minutes for the Zenodo webhook.
3. Set up Python.
4. Install `doi2cff`.
5. Convert the latest TARDIS release DOI to `CITATION.cff`.
6. Try 10 times with a 60 second sleep between attempts.
7. Upload `CITATION.cff` as an artifact.

Credits job:

1. Check out TARDIS.
2. Wait 3 minutes for the Zenodo webhook.
3. Set up Python.
4. Install `requests`.
5. Run a helper script to update `README.rst` and `docs/resources/credits.rst`.
6. Upload `README.rst` and `credits.rst` as artifacts.
7. Dispatch updates to the TARDIS website.

Post-release pull request job:

1. Check out TARDIS.
2. Download artifacts from previous steps.
3. Copy `CHANGELOG.md`, `CITATION.cff`, `README.rst`, and `credits.rst` into
   the repository.
4. Get the current date.
5. Create a pull request.
6. Wait 30 seconds for the pull request to be created.
7. Automatically approve the pull request using tokens from infrastructure and
   core coordinator members.
8. Enable auto-merge.

##### TARDIS Carsus Compatibility Check

The TARDIS Carsus Compatibility Check, or Bridge, compares reference data
generated with different Carsus versions. It has two jobs:

- `carsus-build`: generates an atomic file with the latest Carsus.
- `tardis-build`: generates new reference data with the atomic file.

The two reference data files are compared with this notebook:

https://github.com/tardis-sn/tardis-refdata/blob/master/notebooks/ref_data_compare_from_paths.ipynb

The workflow has a `workflow_dispatch` event for manual runs and is also
triggered weekly by the save atomic files workflow.

##### Save Atomic Files Workflow

The Save Atomic Files workflow runs weekly and can be triggered manually. It
runs the Bridge and sends an artifact to Moria containing the generated atomic
data file and comparison notebook. A separate job indicates whether the Bridge
failed.

##### Regression Data Comparison Workflow

The Regression Data Comparison workflow compares regression data between the
current branch and the base branch on pull requests. It runs only on pull
requests, not on `master`.

The workflow generates regression data for the latest commit on the pull request
and compares it with `master` using the comparison notebook. The notebook is
uploaded as an artifact and pushed to the `reg-data-comp` repository for
previews in the bot comment.

The workflow exports images from the comparison notebook and embeds them in the
bot comment. Unless there are key changes to HDF files in the regression data,
the bot shows two images: spectrum change and relative changes in keys. If
there are key changes, the bot shows a third image visualizing key changes.

##### LFS-Cache Workflow

The `LFS-cache` workflow caches regression data and atomic data. It can be
triggered manually or when there is a push to the main branch of the regression
data repository. It performs LFS pulls when necessary and caches objects. The
`setup-lfs` action restores cached objects. Both fail if the cache is
unavailable.

### Developer Workflow

#### Tutorial: Add a New Feature

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

<span style="color:red">Added: concrete feature-development example.</span>

For example, if a change adds support for a new density configuration, a small
development path would be: open an issue describing the new density option, add
a failing test near the existing model or configuration tests, update
`tardis/io/model/parse_density_configuration.py`, run the targeted tests, and
then document the new configuration behavior.

#### How-To Guide: Introduce a New Feature

Before implementing a new feature, make the feature reviewable:

1. Open or find an issue that states the user need, the expected behavior, and
   any scientific or data assumptions.
2. Identify the smallest code path that owns the behavior. For example, a new
   model input format likely belongs near `tardis/io/model/`, while a new
   visualization behavior likely belongs near `tardis/visualization/`.
3. Add a failing unit test or regression test that captures the new behavior.
4. Implement the smallest code change that makes the test pass.
5. Add or update documentation for the feature.
6. Run targeted tests, Ruff, and any relevant documentation build.
7. Open a pull request with the issue link, test commands, and reviewer notes.

Example feature checklist:

```markdown
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
```

<span style="color:red">Added: explicit feature-introduction workflow.</span>

#### Explanation: Developer Workflow

The developer workflow is intended for anyone contributing to TARDIS
development. Much of the workflow is taken from Astropy's development workflow,
with credit to the Astropy team.

The workflow emphasizes:

- Filing issues for bugs, inconsistencies, missing functionality, and planned
  features.
- Working from a fork.
- Keeping an `upstream` remote connected to the main TARDIS repository.
- Avoiding development on local `master`.
- Starting each separable task from current upstream trunk.
- Keeping work organized through feature branches.
- Running tests before review.
- Updating documentation.
- Opening pull requests for review and merge.

This helps maintain readable history and helps maintainers understand what
changed and why.

#### Explanation: Test-Driven Development

TARDIS aims to use test-driven development. The preferred path starts by
describing a feature, then writing a test for the missing behavior, then
implementing the feature until the test passes.

Unexpected test failures can indicate that the test is wrong. A passing test
before implementation can mean the feature already exists and only lacked
coverage.

### External Links

#### Reference: External Link Reference

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

<span style="color:red">Deleted: removed several older general-purpose Git
learning links from this draft to keep the reference focused on links directly
used by the developer workflow.</span>

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

### Git Workflow

#### How-To Guide: Connect Your Fork to Upstream

Add the main TARDIS repository as `upstream`:

```shell
git remote add upstream https://github.com/tardis-sn/tardis.git
```

The `upstream` name is the convention used for the main TARDIS repository. The
workflow uses the HTTPS URL for upstream because it is read-only by default for
most contributors and avoids accidental writes to the main repository.

Confirm the remotes:

```shell
git remote -v show
```

Expected output should include:

```text
upstream   https://github.com/tardis-sn/tardis.git (fetch)
upstream   https://github.com/tardis-sn/tardis.git (push)
origin     git@github.com:your-user-name/tardis.git (fetch)
origin     git@github.com:your-user-name/tardis.git (push)
```

#### How-To Guide: Fork and Clone TARDIS

You only need to fork once for each package you want to contribute to.

1. Log into GitHub.
2. Go to the TARDIS GitHub home page.
3. Click the fork button.
4. Clone your fork:

   ```shell
   git clone git@github.com:your-user-name/tardis.git
   cd tardis
   ```

5. Check branches:

   ```shell
   git branch -a
   ```

6. Check remotes:

   ```shell
   git remote -v
   ```

At this point, `origin` should point to your GitHub fork.

#### Explanation: Git Branching Strategy

In the developer workflow, `master` is referred to as the trunk. Feature work
should begin from `upstream/master`.

Feature branches make reviews easier because each branch contains a related set
of edits. Contributors should avoid merging trunk or other branches into a
feature branch when possible. If trunk changes need to be incorporated, rebasing
is preferred because it replays feature commits on top of the latest trunk and
keeps history easier to read.

A rebase transforms a history like:

```text
      A---B---C cool-feature
     /
D---E---F---G trunk
```

into:

```text
              A'--B'--C' cool-feature
             /
D---E---F---G trunk
```

Merge conflicts can occur if the feature branch and trunk changed the same
files. Resolve conflicts using the Git rebase documentation and related merge
resolution guidance.

#### Reference: Git Command Reference

Clone a fork:

```shell
git clone git@github.com:your-user-name/tardis.git
```

Inspect branches:

```shell
git branch -a
```

Inspect remotes:

```shell
git remote -v
git remote -v show
```

Add upstream:

```shell
git remote add upstream https://github.com/tardis-sn/tardis.git
```

Fetch upstream:

```shell
git fetch upstream
```

Start a branch from trunk:

```shell
git checkout upstream/master
git checkout -b my-new-feature
```

Push a branch:

```shell
git push origin my-new-feature
git push --set-upstream origin my-new-feature
```

Check status and diff:

```shell
git status
git diff
```

Stage files:

```shell
git add new_file_name
git add modified_file_name
```

Commit:

```shell
git commit -m "A commit message"
```

Rebase:

```shell
git fetch upstream
git checkout cool-feature
git branch tmp cool-feature
git rebase upstream/master
```

Abort a rebase:

```shell
git rebase --abort
```

Force push a rebased branch to your fork:

```shell
git push -f origin cool-feature
```

Reset to a backup branch:

```shell
git reset --hard tmp
```

Inspect reflog:

```shell
git reflog show cool-feature
```

#### Tutorial: Prepare Your Local Development Fork

This path prepares a new contributor to work on TARDIS locally.

1. Set up a Python environment. The original developer workflow recommends
   Anaconda and refers readers to the installation guide.
2. Create a GitHub account if you do not already have one.
3. Configure your GitHub account for write access, including SSH keys.
4. Fork the TARDIS repository from the TARDIS GitHub home page.
5. Clone your fork:

   ```shell
   git clone git@github.com:your-user-name/tardis.git
   cd tardis
   ```

6. Inspect local and remote branches:

   ```shell
   git branch -a
   git remote -v
   ```

7. Add the main TARDIS repository as `upstream`:

   ```shell
   git remote add upstream https://github.com/tardis-sn/tardis.git
   ```

8. Confirm the remote configuration:

   ```shell
   git remote -v show
   ```

   Expected remotes include `origin`, pointing to your fork, and `upstream`,
   pointing to `https://github.com/tardis-sn/tardis.git`.

9. Install TARDIS in development mode:

   ```shell
   pip install -e .
   ```

   This installs TARDIS so imports use your repository clone regardless of your
   working directory. Edits in your clone are available the next time you start a
   Python interpreter and `import tardis`.

#### How-To Guide: Recover from Git Mistakes

If a rebase goes wrong while it is in progress:

```shell
git rebase --abort
```

If you notice the problem after the rebase and made a backup branch:

```shell
git reset --hard tmp
```

If you forgot to make a backup branch, inspect the reflog:

```shell
git reflog show cool-feature
```

Example reflog:

```text
8630830 cool-feature@{0}: commit: BUG: io: close file handles immediately
278dd2a cool-feature@{1}: rebase finished: refs/heads/my-feature-branch onto 11ee694744f2552d
26aa21a cool-feature@{2}: commit: BUG: lib: make seek_gzip_factory not leak gzip obj
```

Reset to the point before the bad rebase:

```shell
git reset --hard cool-feature@{2}
```

#### How-To Guide: Set Up GitHub Access

1. Create a GitHub account at https://github.com if needed.
2. Configure write access with SSH keys.
3. Use GitHub Help's SSH key instructions:
   https://help.github.com/articles/generating-ssh-keys

#### How-To Guide: Start a Feature Branch

In the developer workflow, the TARDIS `master` branch is the trunk.

1. Do not use your local `master` branch for development. Consider deleting it
   to reduce confusion.
2. Fetch the latest upstream changes:

   ```shell
   git fetch upstream
   ```

3. Start from the current upstream trunk:

   ```shell
   git checkout upstream/master
   ```

4. Create a new feature branch:

   ```shell
   git checkout -b my-new-feature
   ```

Use a new branch for each separable set of changes: one task, one branch. Choose
an informative name, such as `bugfix-for-issue-14`,
`refactor-density-parser`, or `update-regression-data-docs`.

<span style="color:red">Added: replaced generic branch-name examples with
TARDIS-oriented examples.</span>

Push the branch to your fork:

```shell
git push origin my-new-feature
```

With Git 1.7 or newer, you can set the upstream branch:

```shell
git push --set-upstream origin my-new-feature
```

#### How-To Guide: Use the Editing Workflow

1. Make your changes.
2. Run tests to check for regressions:

   ```shell
   pytest tardis
   ```

3. If Sphinx is installed, check the documentation build:

   ```shell
   cd docs
   make html
   ```

   The build should succeed and should not report warnings.

4. Check changed files:

   ```shell
   git status
   ```

5. Inspect the actual changes:

   ```shell
   git diff
   ```

6. Add new files:

   ```shell
   git add new_file_name
   ```

7. Add modified files you want to commit:

   ```shell
   git add modified_file_name
   ```

8. Check what will be committed:

   ```shell
   git status
   ```

9. Commit:

   ```shell
   git commit -m "A commit message"
   ```

10. Push to your fork:

    ```shell
    git push
    ```

<span style="color:red">Added: concrete editing workflow example.</span>

For a small parser change, a typical local loop might be:

```shell
ruff check tardis/io/model/parse_density_configuration.py
pytest tardis/io/model/readers/tests
git diff
git add tardis/io/model/parse_density_configuration.py
git commit -m "Improve density parser validation"
git push
```

### Issues

#### How-To Guide: Report Issues

TARDIS is under constant development, so bugs, inconsistencies, and missing
functionality can occur. File issues on the official GitHub repository:

https://github.com/tardis-sn/tardis

Interested users are encouraged to contribute to TARDIS development by following
the development workflow.

### Matterbridge

#### How-To Guide: Edit the Cron Job

SSH to OpenSupernova.org and run:

```shell
sudo crontab -e
```

#### Explanation: Matterbridge

Matterbridge connects messaging channels across platforms such as Slack,
Mattermost, and Gitter. The `matterbridge` binary is used with a TOML
configuration file:

```shell
./matterbridge -conf config-tardis-matterbridge.toml
```

The TOML file contains parameters required to connect rooms, including tokens
and passwords. When the application runs, messages can be shared between the
connected rooms.

TARDIS keeps a service running on the OpenSupernova.org server to run
Matterbridge as a daemon. Configuration files are stored in a private GitHub
repository, including the custom TOML and Linux service files.

The server runs Ubuntu 14.04, so TARDIS uses an Upstart script instead of a
Systemd service file. A Systemd file is also included in the repository. A cron
job restarts the service periodically to prevent disconnections.

#### Reference: Matterbridge Reference

Matterbridge command:

```shell
./matterbridge -conf config-tardis-matterbridge.toml
```

Matterbridge configuration:

- Use the TOML file in `tardis-matterbridge` as an example.
- Configure as many gateways as needed.
- For Slack, follow the Matterbridge wiki Slack setup steps and read the linked
  comment warning not to add scopes manually.
- For Gitter, create a dedicated GitHub bot account and copy the token.
- For Mattermost, follow the Matterbridge wiki Mattermost setup steps.

Important paths:

- Matterbridge executable: `/usr/local/bin/matterbridge`
- Matterbridge TOML config:
  `/usr/local/etc/matterbridge/config-tardis-matterbridge.toml`
- Upstart service config: `/etc/init/matterbridge.conf`

Important repositories and services:

- OpenSupernova.org server: http://opensupernova.org
- Private Matterbridge repository:
  https://github.com/tardis-sn/tardis-matterbridge
- Matterbridge releases: https://github.com/42wim/matterbridge/releases/latest

#### How-To Guide: Set Up Matterbridge

<span style="color:red">Deleted: removed repeated Matterbridge overview because
the explanation above already defines the service.</span>

1. SSH to the OpenSupernova.org server.
2. Download the Matterbridge binary for Linux from the releases section.
3. Make the file executable and rename it to `matterbridge`.
4. Copy the `matterbridge` executable to `/usr/local/bin`.
5. Clone the `tardis-matterbridge` repository in `$HOME`.
6. Copy `config-tardis-matterbridge.toml` to
   `/usr/local/etc/matterbridge/`.
7. Copy `matterbridge.conf` to `/etc/init/`.
8. Start the service:

   ```shell
   sudo service matterbridge start
   ```

9. Test the gateways.

#### How-To Guide: Update Server Configuration

After updating the TOML file:

1. SSH to the OpenSupernova.org server.
2. Copy the new `config-tardis-matterbridge.toml` to
   `/usr/local/etc/matterbridge/`.
3. Restart the service:

   ```shell
   sudo service matterbridge restart
   ```

4. Test the gateways.
5. If everything works, open a pull request to `tardis-matterbridge` with the
   new TOML file.

### Pull Requests

#### How-To Guide: Open a Pull Request

When you are ready for review or merge consideration:

1. Go to your fork, for example `https://github.com/your-user-name/tardis`.
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

<span style="color:red">Added: example pull request description template.</span>

Example pull request description:

```markdown
Summary
- Update density parsing validation.
- Add a regression or unit test for the new behavior.

Checks
- pytest tardis/io/model/readers/tests
- ruff check tardis/io/model

Notes for reviewers
- Please check whether the error message is clear for invalid density types.
```

#### How-To Guide: Review Pull Requests

GitHub supports line comments and discussion threads on pull requests. Sometimes
it is easier to demonstrate a suggested change by making a pull request on a
pull request.

Add pull request refs to the `upstream` remote in `.git/config`:

```ini
[remote "upstream"]
    url = git@github.com:tardis-sn/tardis.git
    fetch = +refs/heads/*:refs/remotes/upstream/*
    fetch = +refs/pull/*/head:refs/remotes/upstream/pr/*
```

Fetch upstream:

```shell
git fetch upstream
```

You can then check out a pull request branch, for example pull request 116:

```shell
git checkout upstream/pr/116
```

This leaves you in a detached `HEAD` state. Create a branch for your changes:

```shell
git checkout -b helping-with-PR116
```

After making and committing your changes, push them to your fork:

```shell
git push origin helping-with-PR116
```

#### How-To Guide: Update a Pull Request

Pull requests often remain open while feedback is added and other pull requests
are merged. If changes on trunk affect your work, rebase your feature branch on
top of the current trunk.

Fetch trunk and switch to your branch:

```shell
git fetch upstream
git checkout cool-feature
```

Create a temporary backup branch:

```shell
git branch tmp cool-feature
```

Rebase onto upstream trunk:

```shell
git rebase --onto upstream/master upstream/master cool-feature
```

If you are already on `cool-feature`, the command can be shortened:

```shell
git rebase upstream/master
```

When everything looks good, delete the backup branch:

```shell
git branch -D tmp
```

If the feature branch is already on GitHub, force push after rebasing:

```shell
git push -f origin cool-feature
```

Force pushing overwrites the branch on GitHub and can lose commits if used
incorrectly. Never force push to the main TARDIS repository, typically called
`upstream`, because that rewrites shared history.

## Data And Codebase Structure

Guidance for local development setup, codebase-specific design notes, and debugging support.

### Developer FAQ

#### Reference: Developer FAQ

Constants in TARDIS are taken from Astropy. The `tardis.constants` module
imports all constants from `astropy.constants.astropy13constants`.

Class design and inheritance:

- If only the constructor changed, use a classmethod.
- If overriding other methods, use a subclass.

TARDIS uses Ruff to check PEP 8 compliance.

<span style="color:red">Added: examples for classmethod versus subclass
guidance.</span>

Use a classmethod when construction changes but behavior stays the same, such as
adding `from_config(...)` or `from_csvy(...)` constructors around an existing
class. Use a subclass when method behavior changes after construction, such as a
specialized solver or reader that overrides parsing or calculation methods.

### Development Install

#### How-To Guide: Install TARDIS in Development Mode

From the repository root, run:

```shell
pip install -e .
```

TARDIS is designed to be usable from the source tree. An editable install makes
`tardis` import from your clone, so local edits are available immediately in new
Python sessions.

<span style="color:red">Added: import-check example after editable install.</span>

After installation, confirm that Python imports from your checkout:

```shell
python -c "import tardis; print(tardis.__file__)"
```

### Numba Debugging

#### How-To Guide: Debug `numba_montecarlo`

TARDIS provides experimental debugging configurations for deeper debugging when
interfacing with the `montecarlo_numba` module.

PyCharm debugging configurations, related scripts, and `.yml` files are in
`tardis.scripts.debug`. They currently include single-packet mode.

Relevant files:

- `tardis_example_single.yml`: configuration file for the single-packet TARDIS
  run.
- `run_numba_single.py`: Python script that runs the `.yml` file.
- `run_numba_single.xml`: PyCharm debug configuration file used with the files
  above.

This method is experimental.

## Documentation

Guidance for writing, building, previewing, and troubleshooting developer documentation.

### Docstrings

#### Reference: Docstring Reference

A docstring describes a module, function, class, or method definition. It is
stored as `object.__doc__` and is surrounded by triple double quotes.

TARDIS follows the NumPy docstring format:

https://numpydoc.readthedocs.io/en/latest/format.html

Sphinx uses docstrings to auto-generate API documentation:

https://tardis-sn.github.io/tardis/api/modules.html

<span style="color:red">Added: replaced generic docstring sample with a TARDIS
function example.</span>

Example adapted from `tardis/io/model/parse_density_configuration.py`:

```python
def parse_density_section_config(
    density_configuration: ConfigurationNameSpace,
    v_middle: u.Quantity,
    time_explosion: u.Quantity,
) -> tuple[u.Quantity, u.Quantity]:
    """
    Parse the density section of the configuration file.

    Parameters
    ----------
    density_configuration : ConfigurationNameSpace
        Density configuration namespace.
    v_middle : astropy.units.Quantity
        Middle of the velocity bins.
    time_explosion : astropy.units.Quantity
        Time of the explosion.

    Returns
    -------
    density_0 : astropy.units.Quantity
        Density at time_0.
    time_0 : astropy.units.Quantity
        Time of the density profile.
    """
```

Docstring conventions:

- Do not include leading or trailing carriage returns.
- Include a carriage return between each segment.
- Start with a summary of the function, class, module, or method.
- The summary should use standard English syntax, start with a capital letter,
  and end with appropriate punctuation.
- The summary should describe purpose, not individual lines or return values.
- Comments on individual lines should be inline comments.
- Variable, module, function, and class names should be written between single
  backticks in prose.
- In the `Returns` section, always state the type, even if the variable name is
  omitted.
- Do not include a `Returns` section when there is no return value.
- Always list the full path for a variable type if it is not a built-in type,
  such as `astropy.units.Quantity`.

Returns section format:

```python
"""
Returns
-------
(`optional variable name` : )type
    (optional descriptor)
"""
```

### Documentation Builds

#### How-To Guide: Build Documentation Locally

Build the documentation:

```shell
cd docs
make html
```

Notes:

- On a fresh local copy, you may need to run `pip install -e .` first.
- Use `DISABLE_NBSPHINX=1 make html` to disable notebook rendering for a faster
  build.
- Use `make html NCORES=<number of cores>` to build in parallel.
- Use `make html NCORES=auto` to use all available device cores.
- Use `make html SPHINXOPTS="<insert sphinx options>"` to pass additional
  Sphinx options.

After the build, open `docs/_build/html/index.html` in your browser. Navigate to
the changed or added page and check that it looks as intended. Check the
terminal for warnings, often caused by faulty hyperlinks or missing toctree
entries. Fix warnings before merging.

<span style="color:red">Added: example fast documentation build command.</span>

When checking only RST changes or link structure, use:

```shell
cd docs
DISABLE_NBSPHINX=1 make html NCORES=auto
```

#### Reference: Documentation Command Reference

Build documentation:

```shell
cd docs
make html
```

Build without notebook rendering:

```shell
DISABLE_NBSPHINX=1 make html
```

Build in parallel:

```shell
make html NCORES=<number of cores>
make html NCORES=auto
```

Pass Sphinx options:

```shell
make html SPHINXOPTS="<insert sphinx options>"
```

Built docs location:

```text
docs/_build/html/index.html
```

Pull request documentation preview URL:

```text
https://tardis-sn.github.io/tardis/pull/<pull request number>/index.html
```

Notebook file-size workaround:

```python
%config InlineBackend.figure_formats='png2x'
```

#### How-To Guide: Share Built Documentation in a Pull Request

Add the `build-docs` label to your pull request. If you cannot add the label,
leave a comment in the pull request or contact a senior member of the
collaboration.

The documentation builds when the label is added. Subsequent commits trigger new
documentation builds while the label remains present.

The built documentation is available at:

```text
https://tardis-sn.github.io/tardis/pull/<pull request number>/index.html
```

It is also linked automatically in pull request comments.

To view build logs, go to the Actions tab in the TARDIS repository and select
the `docs` workflow. Search documentation builds by branch to find the relevant
log.

#### How-To Guide: Troubleshoot Documentation Builds

Documentation should be free of warnings and errors.

Common causes:

- Errors often mean notebooks are incompatible with new code. Update notebooks
  to reflect the code changes.
- Warnings often come from incorrect RST syntax in links, section headers,
  tables of contents, or similar structures.
- Warnings can also come from docstrings that do not follow the NumPy docstring
  format.
- GitHub built documentation files, including `.ipynb` files built by Sphinx,
  can be at most 100 MB.
- Check file sizes after a local documentation build in `docs/_build/html`.
- Notebook image output built by Sphinx defaults to SVG. Detailed SVG images can
  be very large.
- If file size becomes a problem, add
  `%config InlineBackend.figure_formats='png2x'` in a hidden cell at the
  beginning of the notebook.

The Sublime and Sphinx Guide is a useful resource for RST syntax:

https://sublime-and-sphinx-guide.readthedocs.io/en/latest/index.html

Reach out for help if documentation issues are difficult to resolve.

### Documentation Pages

#### Reference: Developer Documentation Map

General developer workflow pages:

- Reporting issues
- Git workflow
- Documentation guidelines
- Running tests
- Benchmarks
- Code quality
- Developer FAQ

Advanced core team development pages:

- Continuous integration
- Updating regression data
- Matterbridge
- Debugging `numba_montecarlo`

#### How-To Guide: Document Code Changes

When you make or add functionality in TARDIS, create or update an `.rst` file or
Jupyter notebook (`.ipynb`) to demonstrate how the feature works. Include the
page in the documentation.

For RST documentation:

1. Use reStructuredText for pages that do not feature interactive code examples.
2. Commit the `.rst` source file, not built HTML.
3. Keep the documentation clear and concise.

For notebook documentation:

1. Use Jupyter notebooks when code examples help explain concepts.
2. TARDIS uses `nbsphinx` to turn notebooks into HTML pages.
3. During documentation builds, `nbsphinx` runs notebooks with cleared output and
   places generated output in the HTML.
4. Always clear notebook output before submitting notebooks.
5. In VS Code, use the "Clear All Outputs" command.
6. In JupyterLab, use `Edit > Clear Outputs of All Cells`.

Running notebooks during the documentation build helps keep documentation
up to date. If code updates are inconsistent with documentation, the
documentation build returns an error.

TARDIS notebook documentation can provide interactive tutorials. Documentation
notebooks include a launch button that directs readers to Binder, where the
notebook can run with an online Jupyter kernel. Notebooks in the Input/Output
section are automatically linked from the tutorials page.

When new functions or classes are added, add docstrings as well as page-level
documentation. Sphinx uses docstrings to auto-generate the API documentation for
the TARDIS package. Build the documentation and check how the corresponding
module API looks.

#### Explanation: Documentation System

High-quality and consistent documentation helps users find specific tasks and
helps developers understand best practices.

TARDIS uses Sphinx to generate documentation. Sphinx translates source files,
often reStructuredText or Jupyter notebooks, into HTML files and automatically
produces cross-references and indices. Developers new to Sphinx should read the
Sphinx quickstart guide.

RST documentation is used for pages without interactive code examples. Notebook
documentation is used when code examples help explain concepts. TARDIS uses
`nbsphinx` to render notebooks into HTML.

Only source documentation files are committed. Built HTML is not committed.

Notebook output must always be cleared before submission because `nbsphinx` runs
notebooks during the build and inserts fresh output into the generated HTML.
This helps keep documentation synchronized with the current codebase.

#### How-To Guide: Include a New Documentation Page

Whether the page is reStructuredText or a Jupyter notebook:

1. Determine the appropriate location in the documentation. Ask the TARDIS
   collaboration for help if needed.
2. Place the file in the corresponding directory under `docs/`.
3. Include the file in a `toctree` in the corresponding `index.rst`.

Example: a page under "Visualization Tools and Widgets" in the Input/Output
section belongs in a corresponding directory under `docs/` and must be included
in that section's `index.rst`.

<span style="color:red">Added: concrete toctree example.</span>

For example, a new developer page at
`docs/contributing/development/new_tooling.rst` would be linked from
`docs/contributing/development/index.rst` like this:

```rst
.. toctree::
    :maxdepth: 2

    issues
    git_workflow
    documentation_guidelines
    new_tooling
```

## Testing And Validation

Guidance for tests, regression data, benchmarks, and validation workflows.

### Benchmarks

#### Explanation: Benchmarking

The benchmarking system detects performance regressions in TARDIS. It lets
developers visually check whether performance has improved or worsened.

TARDIS uses AirSpeed Velocity, or ASV. ASV is designed to run benchmarks on
random servers, such as GitHub-hosted runners, and reduce noise caused by
technical differences between servers. ASV produces graphs that indicate
whether a regression occurred and can identify commits that affected performance
in specific functions.

Benchmark files live in the `benchmarks/` directory. Results are stored under
`.asv/`.

#### Reference: Benchmark Command Reference

<span style="color:red">Deleted: condensed repeated ASV setup commands because
the benchmark how-to already gives the full sequence.</span>

Common ASV commands:

```shell
asv setup
asv machine --yes
asv run
asv publish
asv preview
```

<span style="color:red">Added: example benchmark location in TARDIS.</span>

Benchmark classes live under `benchmarks/`; for example,
`benchmarks/spectrum_formal_integral.py` defines
`BenchmarkTransportMontecarloFormalIntegral` with ASV methods such as
`time_intensity_black_body`.

#### How-To Guide: Run Benchmarks

TARDIS uses AirSpeed Velocity, or ASV, for benchmarks.

Install ASV and its required environment tooling. ASV needs Conda or Miniconda
and Mamba. Mini-forge includes these installers and can simplify configuration.

Create the benchmark environment:

```shell
export MAMBA_ENV_NAME="benchmark"
mamba create --yes --name "${MAMBA_ENV_NAME}" python asv mamba
mamba init
```

Set up ASV for TARDIS:

```shell
cd tardis
export MAMBA_ENV_NAME="benchmark"
mamba activate "${MAMBA_ENV_NAME}"
asv setup
asv machine --yes
```

Run and publish benchmarks:

```shell
cd tardis
export MAMBA_ENV_NAME="benchmark"
mamba activate "${MAMBA_ENV_NAME}"
asv run
asv publish
```

Preview benchmark output:

```shell
asv preview
```

You can also view the generated data with a local web server of your choice.

<span style="color:red">Added: example single-benchmark workflow.</span>

When iterating on one benchmark file, inspect the benchmark name in
`benchmarks/spectrum_formal_integral.py`, then run a matching ASV benchmark
locally before publishing results:

```shell
asv run -b time_intensity_black_body
asv preview
```

### Regression Data

#### Explanation: Regression Data

Regression data tests compare TARDIS output against stored reference data. The
data lives in the `tardis-regression-data` repository and may include formats
such as NPY, HDF5, and CSV.

The project uses Git LFS for regression data. The regression data repository
should be cloned outside the main `tardis` repository.

The original documentation notes that TARDIS migrated from `tardis-refdata` to
`tardis-regression-data`.

#### Reference: Regression Data Push Error

If pushing regression data fails with an error similar to:

```text
remote: error: GH008: Your push referenced at least 1 unknown Git LFS object:
remote:     <file-hash>
remote: Try to push them with 'git lfs push --all'.
```

Check the LFS endpoint:

```shell
git lfs env
```

Expected endpoint:

```text
Endpoint=https://tardis:tardis-2025-lfs@registry.moria.egr.msu.edu/repository/tardis-lfs/info/lfs (auth=basic)
```

If the endpoint is different, check `.git/config` and `.lfsconfig`. Git
prioritizes `.git/config` over `.lfsconfig`, so remove duplicates and make sure
`.lfsconfig` matches the TARDIS Regression Data repository. See the
`git-lfs-config` documentation:

https://github.com/git-lfs/git-lfs/blob/main/docs/man/git-lfs-config.adoc

#### How-To Guide: Run Regression Data Tests

Advanced tests require TARDIS regression data from:

https://github.com/tardis-sn/tardis-regression-data

<span style="color:red">Deleted: removed repeated migration and repository
format explanation already covered in the regression-data explanation.</span>

Clone the regression data outside the main `tardis` repository:

```shell
GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/tardis-sn/tardis-regression-data.git
cd tardis-regression-data
git lfs fetch
git lfs checkout
```

`GIT_LFS_SKIP_SMUDGE=1` defers downloading large files until the explicit
`git lfs fetch` command. This makes download progress clearer and shows which
files are being retrieved.

Run regression tests:

```shell
pytest tardis --tardis-regression-data=/path/to/tardis-regression-data/
```

Run tests for a specific file or directory:

```shell
pytest tardis/path/to/test_file_or_directory --tardis-regression-data=/path/to/tardis-regression-data/
```

<span style="color:red">Added: concrete regression-test path example.</span>

For example, to run model-reader tests with local regression data:

```shell
pytest tardis/io/model/readers/tests --tardis-regression-data=/path/to/tardis-regression-data/
```

Some cases require updating the regression data. See
[Update Regression Data](#how-to-guide-update-regression-data).

The tests workflow runs on pull requests and push events. To prevent leaking LFS
quota, tests are disabled on forks. If you need to run tests on a fork, run the
tests workflow on the `master` branch first. The LFS cache generated on `master`
should be available in child branches. Check whether the cache was generated in
the `Setup LFS` step of the workflow run. The cache can also be found under the
Management section of the Actions tab.

#### How-To Guide: Update Regression Data

Regression data tests run only when `pytest` is called with the
`--tardis-regression-data` flag. These tests compare TARDIS output, mostly
arrays, against information stored in regression data files.

If regression data tests fail in a pull request, possible causes include:

A. There is a problem in your code.  
B. Your code is correct, but the regression data is outdated.  
C. The pipeline is broken.

If you suspect the regression data is outdated:

1. Activate the `tardis` environment.
2. Fork and clone the `tardis-regression-data` repository.
3. Follow any instructions inside your local copy of `tardis-regression-data`.
4. Go to your local `tardis` repository.
5. Make sure you are on the branch from which you want to generate new
   regression data.
6. Generate new regression data:

   ```shell
   pytest tardis --tardis-regression-data=/path/to/tardis-regression-data --generate-reference
   ```

7. Check the results and make sure everything is correct.
8. Make a new branch in `tardis-regression-data`.
9. Push the new regression data.
10. Open a pull request.

If issues arise, tag a TARDIS team member responsible for CI/CD:

https://tardis-sn.github.io/people/collaboration/

<span style="color:red">Added: example from the regression comparison workflow.</span>

The pull-request regression comparison workflow uses the same idea in
`.github/workflows/compare-regdata.yml`:

```shell
pytest tardis ${PYTEST_FLAGS} --generate-reference -m "not continuum"
```

### Test Coverage

#### How-To Guide: Generate and View Test Coverage

TARDIS already uses coverage in CI. The test workflows pass coverage flags such
as `--cov=tardis`, `--cov-report=xml`, and `--cov-report=html`, then combine
coverage artifacts with `coverage combine`, `coverage xml`, and `coverage html`.
The older `tox.ini` coverage environment also writes `coverage.xml`.

To generate local coverage for a targeted area:

```shell
pytest tardis/io/model/readers/tests --cov=tardis --cov-report=xml --cov-report=html
```

Then inspect the generated report:

```shell
python -m http.server --directory htmlcov 8000
```

Open `http://localhost:8000` in a browser to inspect line-by-line coverage.

For VS Code, a popular option is the Coverage Gutters extension
(`ryanluker.vscode-coverage-gutters`). It displays coverage from XML or LCOV
files in the editor gutter.

Install it from VS Code Quick Open:

```text
ext install ryanluker.vscode-coverage-gutters
```

Basic workflow:

1. Generate `coverage.xml` with pytest coverage.
2. Open the TARDIS repository in VS Code.
3. Run `Coverage Gutters: Watch` from the command palette.
4. Open a Python file. Covered, uncovered, and partially covered lines appear
   in the gutter.
5. Re-run tests and refresh coverage when the report changes.

See:

- Coverage Gutters marketplace page:
  https://marketplace.visualstudio.com/items?itemName=ryanluker.vscode-coverage-gutters
- Coverage.py command reference:
  https://coverage.readthedocs.io/en/latest/cmd.html

<span style="color:red">Added: test coverage workflow and VS Code coverage
viewer instructions.</span>

### Tests

#### How-To Guide: Run Unit Tests

TARDIS focuses primarily on unit tests. These test individual functions and run
quickly.

Run the unit tests from the repository root:

```shell
pytest tardis
```

<span style="color:red">Added: targeted pytest examples.</span>

Run one package's tests:

```shell
pytest tardis/io/model/readers/tests
```

Run one test file:

```shell
pytest tardis/tests/test_util.py
```

#### Explanation: Good and Bad Test Cases

A good test case is narrow, deterministic, and tied to a behavior that matters.
It should fail for the bug or missing feature it describes, and it should not
depend on unrelated runtime state.

Good test case patterns:

- Tests one behavior at a time.
- Uses the smallest realistic fixture or input.
- Has a clear assertion.
- Uses `numpy.testing` helpers for numerical arrays.
- Names the expected behavior in the test name.

Example:

```python
import numpy.testing as npt


def test_ascii_density_reader_converts_velocity_units():
    time_model, velocity, _ = read_simple_ascii_density(...)

    npt.assert_allclose(velocity[3].value, 1.3e4 * 1e5)
    npt.assert_almost_equal(time_model.to(u.day).value, 1.0)
```

A bad test case is broad, fragile, or unclear. It may pass while the behavior is
wrong, or fail because of unrelated changes.

Poor test case patterns:

- Tests many unrelated features at once.
- Uses a full simulation when a small parser/unit test would work.
- Asserts only that code ran without checking the result.
- Uses exact equality for floating-point arrays.
- Depends on local paths, network access, or test ordering.

Poor example:

```python
def test_model():
    result = run_everything()
    assert result is not None
```

<span style="color:red">Added: good-versus-bad test case guidance and examples.</span>

#### Reference: Testing Command Reference

Run unit tests:

```shell
pytest tardis
```

Clone regression data:

```shell
GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/tardis-sn/tardis-regression-data.git
cd tardis-regression-data
git lfs fetch
git lfs checkout
```

Run regression tests:

```shell
pytest tardis --tardis-regression-data=/path/to/tardis-regression-data/
```

Run regression tests for a path:

```shell
pytest tardis/path/to/test_file_or_directory --tardis-regression-data=/path/to/tardis-regression-data/
```

Generate reference data:

```shell
pytest tardis --tardis-regression-data=/path/to/tardis-regression-data --generate-reference
```

#### Explanation: Testing Strategy

TARDIS primarily uses unit tests. Unit tests check individual functions and run
quickly enough to provide immediate feedback after changes.

The tests use `pytest` and are based on the `astropy-setup-helpers` package.

Unit tests are run after suggested changes to TARDIS to maintain code quality
and prevent regressions.
