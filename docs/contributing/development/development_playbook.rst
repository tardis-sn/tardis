.. _development-playbook:

********************************
TARDIS Development Playbook
********************************

This page documents the day-to-day developer workflow used by the TARDIS team.
It complements (and does not replace) the existing installation, testing,
and CI pages.

Source of Truth for Setup
=========================

For environment creation and first-time installation commands, follow
:ref:`installation`.

.. note::

    The installation guide is the standard source for exact bootstrap commands.
    Keep setup instructions centralized there, and link from development pages
    to avoid drift.

Developer QA Loop Before Review
================================

For a typical code pull request, contributors are expected to:

1. Run Ruff formatting/linting (through VS Code integration or CLI):

   .. code-block:: bash

      ruff format .
      ruff check .

2. Run tests locally where possible. For full regression-backed execution:

   .. code-block:: bash

      pytest tardis --tardis-regression-data=/path/to/tardis-regression-data

3. Build documentation locally for documentation-heavy or user-facing changes:

   .. code-block:: bash

      cd docs
      make html

4. Wait for CI completion before merge.

.. note::

    Running tests from VS Code test tooling is encouraged when available.
    CLI commands remain the fallback and should be considered authoritative.

Change-Type Testing Expectations
================================

Use this matrix to decide what to run before opening or updating a PR:

* Documentation-only change: local docs build is preferred; CI docs build can
  be used as fallback.
* Code change: local tests are expected where possible, then CI validation.
* Regression-data-related change: use CI-based regression comparison and open
  matching updates in the regression-data repository when needed.
* CI/workflow change: validate in a fork when feasible, then verify in CI.

Regression Data Policy
======================

General policy:

* Update existing regression data only when behavior changes are intentional,
  understood, and justified by a feature or deprecation path.
* Prefer adding new regression cases instead of replacing existing references.
* Regression data artifacts are versioned in the dedicated
  ``tardis-regression-data`` repository, not in this repository.
* Domain experts and lead maintainers are expected to approve substantive
  regression-data updates.

PRs should be blocked when regression changes are not scientifically acceptable
or not sufficiently justified to reviewers.

CI Triage Guidance
==================

General triage priorities:

* Test failures are required to resolve before merge.
* Documentation and codestyle checks are strongly preferred to pass.
* Coverage and benchmark signals are advisory unless a reviewer explicitly
  elevates them for the PR.
* Known failures caused by upstream dependencies may be treated as non-blocking
  when maintainers confirm they are unrelated to the proposed changes.

Rerun CI only when there is a plausible transient cause (for example,
cache/network flakiness). Otherwise, identify and fix root cause first.

Merge Readiness and Review Expectations
=======================================

Typical merge requirements:

* Clear PR title and description.
* Required checks passing (see below).
* At least one maintainer review, or multiple reviews including less
  experienced maintainers.
* Documentation preview available for documentation-impacting PRs when possible.

Currently emphasized required checks:

* ``tests / tests not continuum ubuntu-latest pip tests disabled``
* ``tests / tests not continuum macos-latest pip tests disabled``
* ``docs / build-docs``

.. note::

    Required-check policy can evolve. Keep this list synchronized with branch
    protection settings and CI workflow ownership decisions.

External Contributors vs Core Maintainers
=========================================

External contributors are expected to follow the full documented workflow,
including passing required checks and providing complete PR context.

Core maintainers may apply limited discretion for urgent fixes (for example,
expedited merges), but should still preserve strict standards for test health
and scientific correctness.

Notebook Expectations
=====================

Adding a new notebook is optional and case-by-case. In practice, notebooks are
most useful when a change introduces user-facing behavior, a new workflow, or a
concept that benefits from executable examples.
