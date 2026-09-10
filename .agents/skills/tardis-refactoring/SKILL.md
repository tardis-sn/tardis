---
name: tardis-refactoring
description: Refactor TARDIS internals while preserving intended scientific and user-observable behavior, without treating Python API compatibility as a constraint. Use for behavior-preserving redesigns, code movement, renaming, signature changes, dependency cleanup, decomposition, consolidation, or removal of obsolete internal interfaces. Do not use alone when scientific results, configuration behavior, file formats, or other user-observable behavior should change; use $tardis-behavior-change for that scope.
---

# Refactor TARDIS

1. Read the root `AGENTS.md` and the current developer documents it routes for
   code style, design, and testing. Read other routed documents only when they
   apply to the affected path. Inspect the owning code, callers, tests,
   fixtures, and documentation before editing.

2. State the design problem, intended structure, and invariants to preserve.
   Distinguish Python interface shape from observable behavior:

   - Preserve scientific results, accepted user inputs, generated outputs, and
     documented behavior unless the request explicitly changes them.
   - Do not preserve functions, classes, methods, imports, signatures, or
     module locations merely for API compatibility. TARDIS has no public API.
   - Update every in-repository production caller, example, and documentation
     reference in the same change. Do not add compatibility aliases, wrappers,
     or deprecation shims unless a concrete in-scope consumer requires them.

3. Establish a passing baseline with the narrowest existing tests that cover
   the invariants. Treat those tests as a fixed behavioral contract and do not
   modify, remove, or replace them by default. Modify a test only when the
   refactor cannot proceed otherwise, such as when it directly imports a moved
   internal interface, constructs an intentionally changed signature, or
   asserts an implementation detail that no longer exists. Before editing it,
   state why a production-only change is impossible. Keep the edit to the
   smallest mechanical migration; do not weaken assertions or change expected
   scientific results.

   Add a focused characterization test only when an important invariant is
   unprotected and preserving it cannot be demonstrated with existing tests.
   Do not add tests merely to facilitate the refactor. A pure refactor uses a
   green-to-green workflow, not an artificial failing test.

4. Make the smallest coherent structural change that solves the stated design
   problem. Remove superseded code and migrate call sites completely. Avoid
   mixing cleanup or behavior changes unrelated to the refactoring goal.

5. If the work reveals that an observable or scientific behavior must change,
   separate that scope and use `$tardis-behavior-change` before implementing
   it. Use `$tardis-regression-validation` when the affected path has
   regression-backed coverage.

6. Run the closest affected tests and Ruff in the `tardis` Conda environment,
   then expand validation according to root `AGENTS.md`. Compare relevant
   numerical or tabular outputs when structural tests alone would not protect
   scientific equivalence. Update and build documentation when references or
   documented interfaces changed.

7. In the handoff, report the design improvement, interfaces intentionally
   changed or removed, preserved behavioral and scientific invariants, files
   changed, any test edits and why each was unavoidable, validation commands
   and outcomes, and every omitted check with its reason.
