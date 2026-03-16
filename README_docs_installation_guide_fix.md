# Installation/Setup Documentation Fix README

## :pencil: Description

**Type:** :beetle: `bugfix`

This change fixes documentation issues reported in setup and installation guides:

- Installation page mismatch and stale URL usage
- Broken quickstart link from installation page
- Developer workflow missing explicit `tardisbase` install guidance
- Developer workflow command that navigated to home directory instead of project directory

### What changed

- In `docs/getting_started/installation.rst`:
  - Added an explicit section anchor for lockfile installation: ``.. _install-with-lockfiles:``
  - Replaced hardcoded external link:
    - `https://tardis-sn.github.io/tardis/installation.html#install-with-lockfiles`
  - With internal Sphinx reference:
    - ``:ref:`here <install-with-lockfiles>` ``
  - Updated quickstart link from:
    - ``Quickstart for TARDIS <quickstart.ipynb>``
  - To:
    - ``:doc:`Quickstart for TARDIS <../quickstart>` ``

- In `docs/contributing/development/git_workflow.rst`:
  - Fixed installation guide path from:
    - ``:doc:`Installation guide <../../installation>` ``
  - To:
    - ``:doc:`Installation guide <../../getting_started/installation>` ``
  - Corrected command from:
    - `cd`
  - To:
    - `cd tardis`
  - Updated develop install command from:
    - `pip install -e .`
  - To:
    - `pip install -e ".[tardisbase,viz]"`

**Fixes:** #3450

## :robot: AI Usage Statement

AI tools were used for this contribution.

- Tool: Codex (GPT-5 coding assistant)
- Usage: Assisted with identifying stale links/commands, applying documentation updates, and drafting this README.
- Verification: I reviewed and understand all AI-assisted changes and can explain them.

## :vertical_traffic_light: Testing

How did you test these changes?

- [ ] Testing pipeline
- [x] Other method (describe)
- [ ] My changes can't be tested (explain why)

### Local verification

- Ran:
  - `rg` checks for stale/broken patterns in updated docs (`installation.html`, old doc paths, `quickstart.ipynb`, stray `cd`, `pip install -e .`)
- Result:
  - No remaining matches for the corrected patterns in edited files.
- Attempted:
  - `python -m sphinx -b dummy -q docs /tmp/tardis-docs-dummy`
- Result:
  - Could not run docs build locally in this environment: `No module named sphinx`.

## :ballot_box_with_check: PR Checklist Mapping

- [x] One pull request for one specific task (docs setup/install fixes only)
- [x] Explicit AI usage statement included
- [x] PR description content prepared in template style
- [ ] PR title includes `[GSoC]` (apply when opening PR)
- [ ] Address review/bot comments within 48 hours (follow during review)
- [ ] Subject understanding confirmed before opening PR
