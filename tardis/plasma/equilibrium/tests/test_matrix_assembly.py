import numpy as np
import numpy.testing as npt
import pandas as pd

from tardis.plasma.equilibrium.matrix_assembly import (
    assemble_bound_bound_rate_matrix,
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


def test_assemble_bound_bound_matrix_adds_direct_fixed_rates() -> None:
    """Fixed rates enter directly without Sobolev scaling."""
    common_arguments = (
        2,
        np.array([0], dtype=np.int64),
        np.array([1], dtype=np.int64),
        np.array([10.0]),
        np.array([2.0]),
        np.array([3.0]),
        np.array([0], dtype=np.int64),
    )

    trapped_matrix = assemble_bound_bound_rate_matrix(
        *common_arguments,
        np.array([0.0]),
    )
    escaping_matrix = assemble_bound_bound_rate_matrix(
        *common_arguments,
        np.array([0.5]),
    )

    npt.assert_array_equal(trapped_matrix, np.array([[-5.0, 0.0], [5.0, 0.0]]))
    npt.assert_array_equal(
        escaping_matrix,
        np.array([[-10.0, 0.0], [10.0, 0.0]]),
    )
