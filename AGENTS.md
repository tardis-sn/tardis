# TARDIS Agent Instructions

## Scope and Priorities

- Repository: TARDIS spectral-synthesis code.
- Audience: research software engineers with graduate-level software-engineering
  skills and postdoctoral-level physics context.
- Prioritize scientific correctness, reproducibility, test health, and then
  maintainability.
- Make minimal, scoped changes. Do not make unrelated refactors while completing
  a task.

## Source of Truth

Follow the current developer documentation, not superseded documentation paths.

- Environment setup: [docs/getting_started/installation.rst](docs/getting_started/installation.rst)
- Feature and test-driven workflow: [docs/current_developers/feature_workflow.rst](docs/current_developers/feature_workflow.rst)
- Testing, coverage, and regression data: [docs/current_developers/testing_and_validation.rst](docs/current_developers/testing_and_validation.rst)
- Code style and type hints: [docs/current_developers/code_style_and_syntax.rst](docs/current_developers/code_style_and_syntax.rst)
- Design decisions and constants: [docs/current_developers/design_philosophy.rst](docs/current_developers/design_philosophy.rst)
- Documentation and docstrings: [docs/current_developers/documentation.rst](docs/current_developers/documentation.rst)
- CI workflows: [docs/current_developers/continuous_integration.rst](docs/current_developers/continuous_integration.rst)
- Benchmarking: [docs/current_developers/benchmarking.rst](docs/current_developers/benchmarking.rst)

## Environment

- TARDIS is supported on macOS and GNU/Linux and uses Conda lockfiles for
  reproducible environments.
- Check if a `tardis` Conda environment exists and ask the user to use it.
- For a new environment, create the `tardis` Conda environment from the
  platform lockfile, activate it, and install the checkout in editable mode.
- Install the `tardisbase` extra before running regression tests and the `viz`
  extra before running visualization tools.
- Use the project's installed dependency versions; do not substitute physical
  constants or add dependencies when existing TARDIS dependencies suffice.
- Run every test, lint, documentation-build, and Python command that imports
  TARDIS in the `tardis` Conda environment. When the shell is not
  activated, prefix the command with `conda run --no-capture-output -n tardis`:

  ```shell
  conda run --no-capture-output -n tardis pytest tardis/path/to/tests
  conda run --no-capture-output -n tardis ruff check tardis/path/to/code
  cd docs && conda run --no-capture-output -n tardis make html NCORES=auto
  ```

- Before reporting a development tool as unavailable or a validation step as
  blocked, check it in the `tardis` environment. Do not infer availability from
  the active shell's `PATH`.

## Feature Workflow

For any new or changed observable or scientific behavior, use the repo skill
`$tardis-behavior-change`.

- Record the science requirement, expected behavior, and data assumptions. If
  the request identifies an issue, include its identifier in the handoff; do
  not create or modify issues without user authorization.
- Preserve the red-green workflow: demonstrate the requested behavior with a
  narrow failing test, confirm the intended failure, and implement the smallest
  change that makes the test pass. If the test already passes, treat it as
  coverage for existing behavior rather than changing code.
- Run targeted tests and Ruff, use `$tardis-regression-validation` when the
  affected path has regression-backed coverage, and document every user-facing
  change.

Break large features into small, quantifiable goals. Do not create or switch
branches unless the task authorizes it; when authorized, branch from current
upstream trunk and never work directly on local `master`.

## Code and Design

- Follow the Ruff configuration in [pyproject.toml](pyproject.toml): permanent
  rules live there and non-permanent rules live in `.ruff.toml`.
- Use descriptive names for scientific and codebase concepts. Functions and
  methods use verb-noun names; properties and variables use noun names.
- Use `_idx` for an integer array position and `_index` for a dictionary-like
  lookup key.
- Add type hints to every new or touched function definition. Use
  `numpy.typing.NDArray` for concrete NumPy array types and `npt.ArrayLike` for
  broad array-like inputs.
- Use NumPy-style docstrings for new or changed public-facing behavior. Include
  parameter and return types, and omit `Returns` when there is no return value.
- Use physical constants from `tardis.constants`, which is based on Astropy's
  `astropy13constants`.
- Introduce Numba only for a numerical hotspot with a stated performance need.
  Verify that the changed code compiles in nopython mode. Leave performance
  measurement to remote CI and report available CI benchmark evidence;
  otherwise use the existing non-Numba implementation.
- Use composition instead of inheritance. Use a classmethod when only
  construction changes and a subclass when method behavior changes.
- Add input validation only for an error that is both likely and directly
  actionable by the user. Let other invalid scientific inputs fail naturally.
- Reuse existing dependencies and pytest fixtures. Do not add single-use helper
  functions, trivial wrapper functions, commented-out code, or test helpers
  when an existing fixture is suitable.

## Testing and Regression Data

- Write tests that are narrow, deterministic, and assert the requested
  behavior.
  Use `numpy.testing` and `pandas.testing` comparison helpers for numerical and
  tabular results.
- Run the closest test module first:

  ```shell
  conda run --no-capture-output -n tardis pytest tardis/path/to/tests
  conda run --no-capture-output -n tardis ruff check tardis/path/to/changed_code
  ```

- Run the broader suite after targeted validation unless it cannot complete in
  the available environment or time budget. State the reason in the handoff if
  it is not run:

  ```shell
  conda run --no-capture-output -n tardis pytest tardis
  ```

- For code exercised by regression-backed tests, use the repo skill
  `$tardis-regression-validation` and run the matching test path with
  `--tardis-regression-data`. Unit-test results are not a substitute.
- Regression data lives in the separate `tardis-regression-data` repository and
  must remain outside this checkout. If the required dataset is unavailable,
  report regression validation as incomplete and state the blocking reason.
- Do not generate, overwrite, or commit regression references unless the task
  explicitly authorizes the behavior change and regression-data update. When
  authorized, inspect generated references, update the paired
  `tardis-regression-data` repository on a separate branch, and document the
  baseline rationale.
- Do not run ASV benchmarks locally. Remote CI owns performance benchmark
  execution. Do not claim a performance result without evidence.

## Documentation and CI

- Update `docs/` in the same change as any user-facing behavior. After changing
  documentation, run the following command and resolve every warning before
  handoff:

  ```shell
  cd docs && conda run --no-capture-output -n tardis make html NCORES=auto
  ```
- Clear outputs from every documentation notebook before committing it.
- Put every new GitHub Actions workflow in `.github/workflows`. In workflows
  that need atomic or regression data, use the shared `setup_lfs` action and
  `lfs-cache` workflow; use `setup_env` and set `bash -l {0}` as the run shell.
- Do not modify CI cache keys or use the `allow_lfs_pull` label without a
  justified need.

## Agent Execution and Handoff

- Before editing, read the applicable documentation and inspect nearby code and
  tests. Preserve all pre-existing user changes in a dirty worktree.
- Do not change regression data, dependencies, CI, release files, or scientific
  assumptions unless the request clearly includes that authority.
- In every handoff, report files changed, scientific or behavioral rationale,
  commands run with outcomes, and every omitted check with its reason.
