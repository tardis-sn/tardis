import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal, assert_series_equal
from scipy.integrate import trapezoid

from tardis.iip_plasma.continuum.base import ContinuumProcess
from tardis.iip_plasma.properties.continuum import (
    integrate_array_by_level_groups,
)


def test_normalize_transition_probabilities_matches_groupby_transform():
    """Compare transition normalization with the pandas groupby baseline."""
    index = pd.MultiIndex.from_tuples(
        [(0, 1), (0, 2), (1, 3), (1, 4), (2, 5)],
        names=["source_level_idx", "destination_level_idx"],
    )
    transition_probabilities = pd.DataFrame(
        {
            "transition_type": [1, 0, 2, -1, 0],
            "lines_idx": [10, 11, -1, 12, 13],
            0: [1.0, 3.0, 0.0, 0.0, 5.0],
            1: [2.0, 2.0, 0.0, 0.0, 15.0],
        },
        index=index,
    )
    expected_probabilities = (
        transition_probabilities.iloc[:, 2:]
        .groupby(level=0)
        .transform(lambda values: values / values.sum())
    )
    expected = pd.concat(
        [transition_probabilities.iloc[:, :2], expected_probabilities],
        axis=1,
    ).fillna(0.0)

    actual = ContinuumProcess._normalize_transition_probabilities(
        transition_probabilities, no_ref_columns=2
    )

    assert_frame_equal(actual, expected)
    assert np.issubdtype(actual["transition_type"].dtype, np.integer)
    assert np.issubdtype(actual["lines_idx"].dtype, np.integer)


def test_integrate_array_by_level_groups_matches_groupby_apply():
    index = pd.MultiIndex.from_tuples(
        [
            (1, 0, 0),
            (1, 0, 0),
            (1, 0, 0),
            (2, 1, 0),
            (2, 1, 0),
            (2, 1, 0),
            (2, 1, 1),
            (2, 1, 1),
        ],
        names=["atomic_number", "ion_number", "level_number"],
    )
    nu = pd.Series([1.0, 1.5, 2.0, 3.0, 3.25, 4.0, 5.0, 6.0], index=index)
    values = np.array(
        [
            [1.0, 2.0, 3.0],
            [1.5, 2.5, 3.5],
            [2.0, 3.0, 4.0],
            [0.5, 1.0, 1.5],
            [0.75, 1.25, 1.75],
            [1.0, 1.5, 2.0],
            [1.25, 1.75, 2.25],
            [1.5, 2.0, 2.5],
        ]
    )

    actual = integrate_array_by_level_groups(values, nu)

    frame = pd.DataFrame(values, index=index)
    frame.insert(0, "nu", nu)
    grouped = frame.groupby(level=[0, 1, 2])
    expected = pd.DataFrame(
        {
            shell: grouped.apply(
                lambda sub, shell=shell: trapezoid(sub[shell], sub["nu"])
            )
            for shell in range(values.shape[1])
        }
    )

    assert_frame_equal(actual, expected)


def test_integrate_array_by_level_groups_matches_series_groupby_apply():
    index = pd.MultiIndex.from_tuples(
        [
            (1, 0, 0),
            (1, 0, 0),
            (1, 0, 0),
            (2, 1, 0),
            (2, 1, 0),
        ],
        names=["atomic_number", "ion_number", "level_number"],
    )
    nu = pd.Series([1.0, 1.5, 2.0, 3.0, 4.0], index=index)
    values = np.array([1.0, 1.5, 2.0, 0.5, 1.0])

    level_groups = tuple(nu.groupby(level=[0, 1, 2]).indices.items())
    actual = integrate_array_by_level_groups(values, nu, level_groups)

    frame = pd.DataFrame({"values": values, "nu": nu}, index=index)
    expected = frame.groupby(level=[0, 1, 2]).apply(
        lambda sub: trapezoid(sub["values"], sub["nu"])
    )

    assert_series_equal(actual, expected)
