"""Verification of read-only diagnostics, not of physical helium rates."""

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis.workflows.nonthermal_helium_diagnostics import (
    canonicalize_underflowed_populations,
    diagnose_level_relaxation,
)


@pytest.mark.parametrize(
    ("parent", "level_values"),
    [
        (0.0, [0.0, 0.0]),
        (1e-21, [0.25e-21, 0.75e-21]),
        (np.ldexp(1.0, -1030), [np.ldexp(1.0, -1032)] * 2),
        (np.nextafter(0.0, 1.0), [np.nextafter(0.0, 1.0), 0.0]),
        (np.finfo(np.float64).tiny, [0.0, 0.0]),
        (1e-310, [0.0, 0.0]),
        (1e-21, [0.0, 0.0]),
        (np.nextafter(0.0, 1.0), [np.nan, 0.0]),
        (np.nextafter(0.0, 1.0), [np.inf, 0.0]),
        (np.nextafter(0.0, 1.0), [-1e-310, 1e-310]),
        (np.nan, [0.0, 0.0]),
        (np.inf, [0.0, 0.0]),
        (-np.nextafter(0.0, 1.0), [0.0, 0.0]),
    ],
)
def test_underflow_policy_preserves_representable_or_invalid_states(
    parent: float, level_values: list[float]
) -> None:
    """Only positive subnormal parents with entirely zero levels may be zeroed.

    Binary-exact nonzero subnormals, normal trace populations, the smallest
    normal boundary, and real contamination are negative controls. Exact
    frame equality (including NaNs) detects accidental abundance floors,
    skipna-based cancellation, or blanket nonfinite replacement.
    """
    level_index = pd.MultiIndex.from_tuples(
        [(28, 8, 0), (28, 8, 1)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    levels = pd.DataFrame(level_values, index=level_index, columns=[36])
    ions = pd.DataFrame(
        [parent],
        index=level_index.droplevel("level_number").unique(),
        columns=[36],
    )
    expected_levels, expected_ions = levels.copy(), ions.copy()
    actual_levels, actual_ions = canonicalize_underflowed_populations(
        levels, ions
    )
    pdt.assert_frame_equal(actual_levels, expected_levels, check_exact=True)
    pdt.assert_frame_equal(actual_ions, expected_ions, check_exact=True)
    assert actual_levels is levels
    assert actual_ions is ions


def test_underflow_policy_repairs_only_represented_zero_groups_by_label() -> (
    None
):
    """Reordered tables preserve ownership and ions without represented levels."""
    level_index = pd.MultiIndex.from_tuples(
        [(28, 8, 1), (2, 0, 0), (28, 8, 0)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    levels = pd.DataFrame(
        [[0.0, 0.0], [1.0, 1.0], [0.0, 0.0]],
        index=level_index,
        columns=[36, 7],
    )
    ions = pd.DataFrame(
        [
            [np.nextafter(0.0, 1.0)] * 2,
            [1.0, 1.0],
            [np.finfo(float).tiny, np.nextafter(0.0, 1.0)],
        ],
        index=pd.MultiIndex.from_tuples(
            [(29, 9), (2, 0), (28, 8)],
            names=["atomic_number", "ion_number"],
        ),
        columns=[7, 36],
    )
    original = ions.copy(deep=True)
    actual_levels, actual_ions = canonicalize_underflowed_populations(
        levels, ions
    )
    expected = original.copy()
    expected.loc[(28, 8), 36] = 0.0
    pdt.assert_frame_equal(actual_ions, expected, check_exact=True)
    pdt.assert_frame_equal(ions, original, check_exact=True)
    assert actual_levels is levels
    # A missing represented parent is not synthesized.
    missing = ions.drop(index=(28, 8))
    _, actual_missing = canonicalize_underflowed_populations(levels, missing)
    pdt.assert_frame_equal(actual_missing, missing, check_exact=True)


@pytest.mark.parametrize("missing_field", ["levels", "parent"])
def test_underflow_policy_preserves_nullable_missing_data(
    missing_field: str,
) -> None:
    """A valid repair must not hide adjacent nullable missing levels or parents.

    Both two-level ions have minimum-subnormal parents. The first has exactly
    zero levels and must be repaired. The second contains a pandas missing
    value: it must remain unchanged, not participate in skip-NA group logic.
    Exact nullable frame equality verifies both preservation and ownership.
    """
    quantum = np.nextafter(0.0, 1.0)
    index = pd.MultiIndex.from_tuples(
        [(28, 8, 0), (28, 8, 1), (29, 9, 0), (29, 9, 1)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    levels = pd.DataFrame(0.0, index=index, columns=[36], dtype="Float64")
    ions = pd.DataFrame(
        quantum,
        index=index.droplevel("level_number").unique(),
        columns=[36],
        dtype="Float64",
    )
    if missing_field == "levels":
        levels.loc[(29, 9), :] = pd.NA
    else:
        ions.loc[(29, 9), 36] = pd.NA
    original_levels, original_ions = levels.copy(), ions.copy()
    expected_ions = ions.copy()
    expected_ions.loc[(28, 8), 36] = 0.0
    actual_levels, actual_ions = canonicalize_underflowed_populations(
        levels, ions
    )
    npt.assert_equal(actual_ions.loc[(28, 8), 36], 0.0)
    if missing_field == "levels":
        npt.assert_equal(actual_ions.loc[(29, 9), 36], quantum)
    else:
        assert pd.isna(actual_ions.loc[(29, 9), 36])
    pdt.assert_frame_equal(actual_levels, original_levels, check_exact=True)
    pdt.assert_frame_equal(actual_ions, expected_ions, check_exact=True)
    pdt.assert_frame_equal(levels, original_levels, check_exact=True)
    pdt.assert_frame_equal(ions, original_ions, check_exact=True)


@pytest.mark.parametrize("level_count", [1, 2, 3, 26, 27])
@pytest.mark.parametrize("at_boundary", [True, False])
def test_underflow_policy_obeys_normalized_allocation_bound(
    level_count: int, at_boundary: bool
) -> None:
    """N zero shares imply parent <= floor(N/2) smallest-subnormal quanta.

    Nearest-even rounding maps each exact share <= half a quantum to zero.
    Summing N nonnegative normalized shares gives the independent parent
    bound. One quantum above that bound cannot all vanish, even for an
    otherwise subnormal parent. Exact comparisons expose off-by-one floors.
    """
    quantum = np.nextafter(0.0, 1.0)
    bound_quanta = level_count // 2
    parent = (bound_quanta + (not at_boundary)) * quantum
    index = pd.MultiIndex.from_tuples(
        [(28, 8, level_number) for level_number in range(level_count)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    levels = pd.DataFrame(0.0, index=index, columns=[36])
    ions = pd.DataFrame(
        [parent], index=index.droplevel("level_number").unique(), columns=[36]
    )
    _, actual = canonicalize_underflowed_populations(levels, ions)
    npt.assert_equal(actual.iloc[0, 0], 0.0 if at_boundary else parent)


@pytest.mark.parametrize("represented_candidate", [False, True])
def test_level_audit_reports_raw_underflow_and_repaired_projection(
    represented_candidate: bool,
) -> None:
    """Distinguish raw candidate and post-reconstruction allocation underflow."""
    level_index = pd.MultiIndex.from_tuples(
        [(28, 8, level_number) for level_number in range(26)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    levels = pd.DataFrame(0.0, index=level_index, columns=[36])
    ions = pd.DataFrame(
        [0.0],
        index=level_index.droplevel("level_number").unique(),
        columns=[36],
    )
    candidate_ions = ions.copy()
    quantum = np.nextafter(0.0, 1.0)
    candidate_ions.iloc[0, 0] = (
        26 * quantum if represented_candidate else quantum
    )
    candidate_levels = levels + (quantum if represented_candidate else 0.0)
    row = diagnose_level_relaxation(
        levels,
        candidate_levels,
        ions,
        candidate_ions,
        ion_relaxation_factors=pd.DataFrame(
            [0.25 if represented_candidate else 0.75], index=[28], columns=[36]
        ),
    ).loc[(28, 8, 36)]
    assert row["classification"] == "underflow_repaired"
    assert bool(row["candidate_underflow_repaired"]) == (
        not represented_candidate
    )
    assert bool(row["relaxed_underflow_repaired"]) == represented_candidate
    npt.assert_equal(
        row["candidate_ion_density"],
        26 * quantum if represented_candidate else quantum,
    )
    npt.assert_equal(row["relaxed_ion_density"], 0.0)
    npt.assert_equal(row["projected_nonfinite_level_count"], 0)


@pytest.mark.parametrize(
    (
        "input_values",
        "candidate_values",
        "input_ion",
        "candidate_ion",
        "classification",
        "projected_bad",
        "introduced_bad",
    ),
    [
        (
            [0.0, 0.0],
            [0.0, 0.0],
            3e-23,
            4e-23,
            "positive_ion_without_level_distribution",
            2,
            2,
        ),
        (
            [np.nan, np.nan],
            [np.nan, np.nan],
            3e-23,
            4e-23,
            "nonfinite_input_and_candidate",
            2,
            0,
        ),
        ([1.0, 0.0], [np.nan, np.nan], 1.0, 1.0, "nonfinite_candidate", 0, 0),
        ([np.nan, np.nan], [1.0, 0.0], 1.0, 1.0, "nonfinite_input", 0, 0),
        ([1.0, 0.0], [0.0, 0.0], 1.0, 0.0, "finite", 0, 0),
        (
            [np.nan, np.nan],
            [np.nan, np.nan],
            0.0,
            0.0,
            "nonfinite_input_and_candidate",
            0,
            0,
        ),
        ([0.0, 0.0], [0.0, 0.0], 0.0, 0.0, "finite", 0, 0),
        ([1.0, 0.0], [1.0, 0.0], 1.0, None, "missing_parent_ion", 2, 2),
        ([1.0, 0.0], [1.0, 0.0], 1.0, np.nan, "nonfinite_parent_density", 2, 2),
    ],
)
def test_level_audit_classifies_empty_and_nonfinite_distributions(
    input_values: list[float],
    candidate_values: list[float],
    input_ion: float,
    candidate_ion: float | None,
    classification: str,
    projected_bad: int,
    introduced_bad: int,
) -> None:
    """Count undefined fractions and distinguish creation from propagation.

    Regime: two Li I levels, shell 53, alpha=0.5. The exact limiting cases
    follow from 0/0, NaN times a positive total, fallback to [1, 0], and the
    zero-ion mask. Counts are integers with no numerical tolerance. A missing
    parent is distinct from an existing zero ion; NaN sums must not be
    described as a physical zero. This is numerical diagnostic verification.
    """
    level_index = pd.MultiIndex.from_tuples(
        [(3, 0, 0), (3, 0, 1)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    ion_index = level_index.droplevel("level_number").unique()
    current_levels = pd.DataFrame(input_values, index=level_index, columns=[53])
    candidate_levels = pd.DataFrame(
        candidate_values, index=level_index, columns=[53]
    )
    current_ions = pd.DataFrame([input_ion], index=ion_index, columns=[53])
    candidate_ions = pd.DataFrame(
        [candidate_ion if candidate_ion is not None else 0.0],
        index=ion_index,
        columns=[53],
    )
    if candidate_ion is None:
        candidate_ions = candidate_ions.iloc[:0]
    snapshots = [current_levels, candidate_levels, current_ions, candidate_ions]
    originals = [frame.copy(deep=True) for frame in snapshots]

    report = diagnose_level_relaxation(*snapshots)

    row = report.loc[(3, 0, 53)]
    assert row["classification"] == classification
    assert row["projected_nonfinite_level_count"] == projected_bad
    assert row["introduced_nonfinite_level_count"] == introduced_bad
    assert bool(row["candidate_parent_present"]) == (candidate_ion is not None)
    if np.isnan(input_values).all():
        assert row["input_nonfinite_level_count"] == 2
        assert np.isnan(row["input_level_sum_strict"])
        npt.assert_equal(row["input_level_sum_skipna"], 0.0)
    for snapshot, original in zip(snapshots, originals, strict=True):
        pdt.assert_frame_equal(snapshot, original)


def test_level_audit_aligns_ion_and_shell_labels() -> None:
    """Locate the undefined Li I distribution despite reversed ion/shell order.

    All levels and parent densities are finite. Only Li I in shell 53 has
    zero levels and a positive parent density, so two new NaNs are expected
    there and nowhere else. Exact integer counts verify label alignment.
    """
    level_index = pd.MultiIndex.from_tuples(
        [(3, 0, 1), (2, 0, 0), (3, 0, 0)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    levels = pd.DataFrame(
        [[0.0, 0.0], [1.0, 1.0], [1.0, 0.0]],
        index=level_index,
        columns=[25, 53],
    )
    ions = pd.DataFrame(
        [[1.0, 1.0], [1.0, 1.0]],
        index=pd.MultiIndex.from_tuples(
            [(3, 0), (2, 0)], names=["atomic_number", "ion_number"]
        ),
        columns=[53, 25],
    )
    report = diagnose_level_relaxation(levels, levels, ions, ions)
    affected = report.loc[report["affected"]]
    assert affected.index.tolist() == [(3, 0, 53)]
    npt.assert_array_equal(affected["introduced_nonfinite_level_count"], [2])


def test_level_audit_uses_actual_adaptive_factor() -> None:
    """A full ion step to zero clears NaNs that half damping would retain."""
    index = pd.MultiIndex.from_tuples(
        [(3, 0, 0), (3, 0, 1)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    levels = pd.DataFrame(np.nan, index=index, columns=[53])
    ions = pd.DataFrame(
        [1.0], index=index.droplevel("level_number").unique(), columns=[53]
    )
    fixed = diagnose_level_relaxation(levels, levels, ions, ions * 0)
    adaptive = diagnose_level_relaxation(
        levels,
        levels,
        ions,
        ions * 0,
        ion_relaxation_factors=pd.DataFrame([1.0], index=[3], columns=[53]),
    )
    npt.assert_equal(
        fixed.loc[(3, 0, 53), "projected_nonfinite_level_count"], 2
    )
    npt.assert_equal(
        adaptive.loc[(3, 0, 53), "projected_nonfinite_level_count"], 0
    )
    npt.assert_equal(adaptive.loc[(3, 0, 53), "relaxed_ion_density"], 0.0)
