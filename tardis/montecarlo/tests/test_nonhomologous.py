import pytest

import tardis.montecarlo.montecarlo_numba.nonhomologous_grid as nonhomologous_grid

@pytest.mark.parametrize(["a", "b", "c", "d", "e", "threshold", "expected_roots"],
        [
        (
            0.0,
            0.0,
            0.0,
            2.0,
            -1.0,
            0.0,
            {"result": [0.5]},
        ),])


def test_quartic_roots(a, b, c, d, e, threshold, expected_roots):
    obtained_roots = nonhomologous_grid.quartic_roots(a, b, c, d, e, threshold)

    assert obtained_roots == expected_roots["result"]

