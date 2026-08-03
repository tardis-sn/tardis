---
name: tardis-behavior-change
description: Implement observable TARDIS behavior changes through the repository's scientific, test-driven workflow. Use for features, behavior-changing bug fixes, changed scientific results, or changed user-facing inputs and outputs that require behavioral tests and documentation. Do not use for read-only analysis, pure documentation edits, or behavior-preserving refactors.
---

# Implement a TARDIS Behavior Change

1. Read `docs/current_developers/feature_workflow.rst` completely and read the
   applicable current developer documentation routed by root `AGENTS.md`.
   Inspect the smallest owning code path, nearby tests, fixtures, and relevant
   user documentation before editing.

2. Record the science requirement, exact expected behavior, scientific or data
   assumptions, and any supplied issue identifier. Split large changes into
   small, independently testable goals. Follow the root authorization boundaries
   for issues, branches, dependencies, CI, releases, regression data, and
   scientific assumptions.

3. Add the narrowest deterministic test that demonstrates the requested
   behavior. Place it near existing coverage, reuse existing fixtures, and use
   NumPy or pandas comparison helpers for numerical or tabular results. Prefer
   explicit inputs, fixtures, and dependency injection over monkey patching.
   Use monkey patching only when a stable external boundary cannot otherwise be
   isolated, and patch that boundary rather than internal implementation
   details.

4. Run the new test in the `tardis` Conda environment and confirm that it fails
   because the behavior is missing. Investigate an unexpected failure before
   changing production code. If the test already passes, retain it as coverage
   and do not change production code merely to force a diff.

5. Implement the smallest change in the code path that owns the behavior.
   Preserve unrelated user changes and apply the scientific and code-design
   rules in root `AGENTS.md`.

6. Re-run the new test, then run the closest affected test module and Ruff on
   the changed paths. Use `$tardis-regression-validation` when the changed path
   has regression-backed coverage; do not substitute unit-test results for
   required regression validation.

7. Update `docs/` for every user-facing behavior change. Build the documentation
   and resolve every warning. Run the broader TARDIS test suite unless the
   available environment or time prevents completion.

8. Do not run ASV benchmarks. Leave performance benchmark execution to remote
   CI and make no performance claim without evidence.

9. Report the following in the handoff:

   - Files changed and the scientific or behavioral rationale.
   - Any supplied issue identifier.
   - Scientific, data, and regression assumptions.
   - Every validation command and its outcome.
   - Documentation-build outcome when applicable.
   - Every omitted or incomplete check and its reason.
