---
name: tardis-tests
description: >
  Design, implement, or review scientific and numerical tests in TARDIS.
  Use whenever a task adds or changes physical behavior, numerical algorithms,
  scientific assertions, regression comparisons, validation evidence,
  tolerances, test fixtures representing physical scenarios, or tests for
  plasma, transport, opacity, spectra, atomic data, Monte Carlo calculations,
  or convergence. Also use when reviewing whether existing tests establish
  scientific correctness. Do not use for purely mechanical formatting,
  packaging, or documentation-only changes.
---

# TARDIS Scientific Test Workflow

## 1. Read context

Read the issue, relevant equations or documentation, the changed production
path, nearby tests, and any subsystem-specific AGENTS.md.

Do not edit files yet.

## 2. Produce a scientific test plan

For each proposed test, state:

- test name;
- scientific claim;
- regime and assumptions;
- minimal physical scenario;
- independent verification;
- why the verification is independent;
- tolerance and justification;
- realistic defect detected;
- evidence class: verification, validation, regression, orchestration,
  or smoke.

Stop and revise any proposal whose only evidence is that a value changed,
is finite, has the expected shape, or matches the implementation under test.

## 3. Choose the verification

Use the strongest practical verification in this order:

1. analytic or exact solution;
2. conservation law or invariant;
3. metamorphic relation;
4. limiting case or convergence behavior;
5. independent implementation or external benchmark;
6. regression reference;
7. structural or smoke assertion.

## 4. Implement readable tests

- Express one scientific claim per test.
- Keep decisive physical inputs and units visible.
- Use scenario-oriented fixtures.
- Make expected equations or relations reviewable in the test body.
- Add Claim, Regime, and Verification to the docstring when not obvious.
- Give assertion failures a scientific interpretation.
- Separate scientific assertions from orchestration assertions.

## 5. Justify tolerances

Trace every tolerance to a solver criterion, conditioning, floating-point
effect, discretization error, convergence order, stochastic uncertainty,
or external uncertainty.

## 6. Verify red and green behavior

Where possible:

1. run the test against the pre-change implementation;
2. verify that it fails because the stated contract is violated;
3. run it against the changed implementation;
4. verify that it passes.

Do not accept a failure caused only by missing APIs, fixture construction,
or unrelated exceptions as the red phase.

## 7. Falsify

Identify at least one realistic defect the test claims to catch. When practical,
apply that defect temporarily, run the narrow test, confirm the intended
assertion fails, and revert the mutation.

## 8. Review independently

Review the completed tests for:

- circular expected values;
- weak change-only assertions;
- unjustified tolerances;
- hidden decisive inputs;
- mixed unrelated claims;
- excessive mocking;
- missing units or index contracts;
- failure to detect plausible scientific defects.

## 9. Report

Produce a claim-to-test table containing:

- test;
- claim;
- regime;
- verification;
- tolerance rationale;
- intended defect;
- red-phase result;
- falsification result;
- evidence class.