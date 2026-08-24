.. _development-playbook:

********************************
TARDIS Development Playbook
********************************

This page documents the day-to-day developer workflow used by the TARDIS team.
It complements (and does not replace) the existing installation, testing,
and CI pages.

The wording on this page is intentionally advisory (recommended and
encouraged) because practices can vary slightly by change type and urgency.

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

2. Run tests locally where possible. For full regression-data testing (described further below on this page):

   .. code-block:: bash

      pytest tardis --tardis-regression-data=/path/to/tardis-regression-data

3. Optional: Build documentation locally (for major documentation changes):

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
* Domain experts and experienced maintainers are expected to approve substantive
  regression-data updates.

Recommended workflow when regression data changes are needed:

1. Open a paired PR in ``tardis-regression-data`` and link it from the code PR.
2. Merge the regression-data PR first.
3. Re-run CI on the code PR and merge only after checks are green.

PRs should be blocked when regression changes are not scientifically acceptable
or not sufficiently justified to reviewers.

.. note:: A specific benefit of validating regression data with our CI pipeline 
  is that we run both a mac and linux architecture github runner. By seeing that 
  both runners pass, we validate the code with both supported architectures, 
  which we have seen to produce small but sometimes impactful numerical differences.

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

If broken or failing CI is ignored, leaving a short comment about why on
the PR is encouraged.

Merge Readiness and Review Expectations
=======================================

Typical merge requirements:

* Clear PR title and description.
* Required checks passing (see below).
* Two approvals, with stale approvals dismissed when new commits are pushed.
* Documentation preview available for documentation-impacting PRs when possible.
* No unresolved review conversations.

Currently emphasized required checks:

* ``tests / tests not continuum ubuntu-latest pip tests disabled``
* ``tests / tests not continuum macos-latest pip tests disabled``
* ``docs / build-docs``

Branch protection and merge behavior
====================================

Current practice:

* Squash merge is used by default.
* Admin overrides and direct pushes are reserved for urgent maintenance
  scenarios.

Examples of urgent scenarios include:

* Unexpected breakage introduced by a recently merged PR.
* Unexpected CI failures affecting tests or docs that require immediate repair.

For complex incidents, creating a retrospective issue is recommended.

Workflow labels used in PRs
===========================

Common labels and effects:

* ``build-docs``: forces docs build/preview when needed.
* ``run-generation-tests``: runs tests focused on regression-data generation.
* ``run-regression-comparison``: compares branch-generated regression data
  against stable regression data.

Suggested usage:

* Apply ``build-docs`` for docs changes not automatically detected by CI, or
  when code changes affect documentation notebooks.
* Apply ``run-generation-tests`` when adding new regression data or when
  generation behavior appears problematic.
* Apply ``run-regression-comparison`` when investigating regression failures or
  preparing likely regression-data updates.

Review cadence and escalation
=============================

We perform a weekly triage of open PRs.

If a PR has no review activity for around two weeks, contributors are
encouraged to post a follow-up comment on the PR to request attention.

External contributors should use GitHub Issues for support requests that do not fit as PR comments.

Commit and changelog notes
==========================

No strict commit-message convention is enforced. Because we squash merge PRs, individual commit messages are lost on master.

The project changelog is generated by CI from ``master`` history, so PR titles
should be written as clear release-note-quality summaries.


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
