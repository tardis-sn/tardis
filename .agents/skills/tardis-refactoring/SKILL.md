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
   - Update every in-repository caller, test, fixture, example, and
     documentation reference in the same change. Do not add compatibility
     aliases, wrappers, or deprecation shims unless a concrete in-scope
     consumer requires them.

3. Establish a passing baseline with the narrowest existing tests that cover
   the invariants. Add a focused characterization test only when important
   behavior is currently unprotected; assert behavior rather than the old
   implementation structure. A pure refactor uses a green-to-green workflow,
   not an artificial failing test.

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
   changed, validation commands and outcomes, and every omitted check with its
   reason.
