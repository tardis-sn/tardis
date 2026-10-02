import numpy as np
import numpy.testing as npt
import pandas as pd

from tardis.plasma.equilibrium.matrix_assembly import (
    construct_rate_matrices,
    sum_rate_frames,
)


def test_sum_rate_frames_reindexes_and_scales() -> None:
    index = pd.MultiIndex.from_tuples(
        [(0, 0), (1, 1)], names=["source", "destination"]
    )
    first = pd.DataFrame([1.0, 2.0], index=index, columns=[0])
    second = pd.DataFrame([3.0], index=index[:1], columns=[0])

    actual = sum_rate_frames([first, second], multipliers=[1.0, 2.0])

    expected = pd.DataFrame([7.0, 2.0], index=index, columns=[0])
    pd.testing.assert_frame_equal(actual, expected)


def test_construct_rate_matrices_sums_duplicate_transitions() -> None:
    index = pd.MultiIndex.from_tuples(
        [(0, 1), (0, 1), (1, 0)], names=["source", "destination"]
    )
    rates = pd.DataFrame([[2.0, 4.0], [3.0, 5.0], [7.0, 11.0]], index=index)

    actual = construct_rate_matrices(rates, (2, 2), "source", "destination")

    npt.assert_allclose(
        actual,
        np.array(
            [
                [[0.0, 7.0], [5.0, 0.0]],
                [[0.0, 11.0], [9.0, 0.0]],
            ]
        ),
    )
