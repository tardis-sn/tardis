from __future__ import annotations

from itertools import pairwise

import numpy as np

from tardis.plasma.properties.continuum_processes.fast_array_util import (
    cumulative_integrate_array_by_blocks,
)


def test_cumulative_integrate_array_by_blocks_matches_python_reference() -> None:
    x = np.array([1.0, 1.5, 2.0, 3.0, 3.25, 4.0, 5.0, 6.0])
    f = np.array(
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
    block_references = np.array([0, 3, 8])

    actual = cumulative_integrate_array_by_blocks(f, x, block_references)

    expected = np.zeros_like(f)
    for shell in range(f.shape[1]):
        for block_start, block_stop in pairwise(block_references):
            cumulative_integral = 0.0
            for row in range(block_start + 1, block_stop):
                cumulative_integral += (
                    (x[row] - x[row - 1])
                    * (f[row, shell] + f[row - 1, shell])
                    / 2.0
                )
                expected[row, shell] = cumulative_integral
            expected[block_start + 1 : block_stop, shell] /= cumulative_integral

    np.testing.assert_allclose(actual, expected, rtol=1e-15, atol=0.0)
