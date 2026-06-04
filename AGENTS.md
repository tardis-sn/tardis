# TARDIS Agent Instructions

## Scope and Audience

- Repository: TARDIS spectral synthesis code.
- Expected audience: research software engineers with graduate-level software engineering skills and postdoctoral-level physics context.
- Priority order for changes: scientific correctness, reproducibility, test health, maintainability.

## Source of Truth

Use existing project docs.

- Environment setup: [docs/getting_started/installation.rst](docs/getting_started/installation.rst)
- Day-to-day QA workflow: [docs/contributing/development/development_playbook.rst](docs/contributing/development/development_playbook.rst)
- Testing details and regression-data setup: [docs/contributing/development/running_tests.rst](docs/contributing/development/running_tests.rst)
- Regression-data update workflow: [docs/contributing/development/update_regression_data.rst](docs/contributing/development/update_regression_data.rst)
- Code quality and docstrings: [docs/contributing/development/code_quality.rst](docs/contributing/development/code_quality.rst)
- Docs authoring/build workflow: [docs/contributing/development/documentation_guidelines.rst](docs/contributing/development/documentation_guidelines.rst)
- CI and required checks: [docs/contributing/development/continuous_integration.rst](docs/contributing/development/continuous_integration.rst)

## Coding Conventions

- Follow Ruff and formatting settings in [pyproject.toml](pyproject.toml).
- Keep physics-facing names readable where scientifically standard; do not rename solely for style if it hurts domain clarity. Prefer descriptive variable names over Greek letters.
- Add or update NumPy-style docstrings for behavior changes.
- Reuse existing fixtures and test patterns from nearby tests before adding new abstractions.
- Use type hinting for all definitions.
- Use the TARDIS constants package for all physical constants. This inherits from a specific version of astropy.

## Testing and Regression Policy

- Write tests before solutions.
- For small logic changes, run targeted tests first, then broader `pytest tardis` where feasible.
- Treat regression-data changes as high scrutiny changes.
- Prefer adding regression cases instead of overwriting existing references unless behavior changes are intentional and justified.
- If a change impacts regression baselines, ensure a paired update in `tardis-regression-data` and document rationale.

## Docs and Notebook Expectations

- New user-facing behavior should usually have docs updates in `docs/`.
- For notebooks committed to docs, clear outputs before commit.
- If docs build warnings/errors appear, fix them before merge unless maintainers explicitly defer.

## Agent Behavior Defaults

- Make minimal, scoped edits; avoid opportunistic refactors unrelated to the task.
- When behavior changes, update tests and docs in the same branch when practical.
- Prefer concrete evidence in summaries: what changed, why, and which tests/docs were run.
- If local tooling differs from docs, follow docs as policy and note local limitations in your final report.
- Do not program defensively. Users of TARDIS are experts and are expected to provide correct inputs.
