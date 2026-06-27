import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

from tardis.iip_plasma.continuum.base import ContinuumProcess


def test_normalize_transition_probabilities_matches_groupby_transform():
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
