import pytest
from numpy.testing import assert_almost_equal

import tardis.transport.montecarlo.nonhomologous_grid as nonhomologous_grid


@pytest.mark.parametrize(
    ["a", "b", "c", "d", "e", "threshold", "expected_roots"],
    [
        (
            0.0,
            0.0,
            0.0,
            2.0,
            -1.0,
            0.0,
            {"result": [0.5]},
        ),
        (
            1.0,
            2.0,
            0.0,
            2.0,
            0.0,
            0.0,
            {"result": []},
        ),
        (1.0, -14.0, 71.0, -154.0, 120.0, 2.5, {"result": [3, 4, 5]}),
    ],
)
def test_quartic_roots(a, b, c, d, e, threshold, expected_roots):
    obtained_roots = nonhomologous_grid.quartic_roots(a, b, c, d, e, threshold)
    obtained_roots.sort()

    assert_almost_equal(obtained_roots, expected_roots["result"])
