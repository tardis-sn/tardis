---
name: tardis-regression-validation
description: Validate TARDIS changes with the external tardis-regression-data repository and safely handle explicitly authorized reference updates. Use when changed code is exercised by tests using regression_data, RegressionData, PlotDataHDF, or tardis_regression_path; when a regression-backed test fails or is skipped; or when the user requests regression-reference validation or updates. Do not use for ordinary unit tests that do not depend on regression data.
---

# Validate TARDIS Regression Behavior

Read `docs/current_developers/testing_and_validation.rst` completely before
acting.

## Prepare the Environment

1. Work from the TARDIS repository root and confirm that the `tardis` Conda
   environment exists.

2. Confirm that the environment provides the regression fixture plugin before
   trusting any pytest result:

   ```shell
   conda run --no-capture-output -n tardis python -c "import tardisbase"
   ```

   If `tardisbase` is missing, install the checkout's `tardisbase` extra only
   when environment modification is authorized. Otherwise report regression
   validation as incomplete. Do not substitute dependency versions.

   ```shell
   conda run --no-capture-output -n tardis \
     python -m pip install -e ".[tardisbase]"
   ```

## Locate the Regression Data

1. Use an existing or user-provided `tardis-regression-data` checkout outside
   the TARDIS checkout. Record its absolute path, commit, and dirty-worktree
   state. Preserve all existing changes.

2. Confirm that required Git LFS objects are materialized. Do not clone the
   regression repository into the TARDIS checkout. If no usable dataset is
   available, report regression validation as incomplete; acquire it only when
   the external download and write are in scope.

## Select and Run Tests

1. Inspect nearby tests and fixture dependencies. Search for direct regression
   use with:

   ```shell
   rg -n "regression_data|RegressionData|PlotDataHDF|tardis_regression_path" \
     tardis/path/to/tests
   ```

   Trace indirect fixtures when the test does not name the regression fixture
   directly.

2. Run the smallest matching test node, file, or directory first:

   ```shell
   conda run --no-capture-output -n tardis pytest \
     tardis/path/to/test_file.py \
     --tardis-regression-data=/path/to/tardis-regression-data \
     -rs
   ```

3. Expand to the affected package and then the broader suite when warranted.
   Treat a skipped regression test as unvalidated, not passing. Unit-test
   results do not replace the matching regression-backed test.

4. Classify a failure as a code defect, an intentionally changed result that
   may need a new baseline, or a dataset or pipeline problem. Do not run
   `--generate-reference` as a diagnostic shortcut.

## Update References Only With Explicit Authorization

1. Generate or overwrite references only when the request explicitly
   authorizes both the behavior change and its regression-data update. Confirm
   that authorization covers using a separate branch in the regression-data
   repository. Inspect and preserve both worktrees and document the scientific
   reason for the new baseline.

2. Generate only the narrow affected path:

   ```shell
   conda run --no-capture-output -n tardis pytest \
     tardis/path/to/test_file.py \
     --tardis-regression-data=/path/to/tardis-regression-data \
     --generate-reference
   ```

3. Inspect every generated change, then rerun the same test without
   `--generate-reference`. Authorization to regenerate does not authorize a
   commit, push, pull request, or unrelated reference update.

## Handoff

Report the matching tests and selection rationale, regression-data path and
commit, commands and outcomes, failures and skips, scientific and data
assumptions, any authorized baseline changes, and every omitted check with its
reason. State explicitly when regression validation remains incomplete.
