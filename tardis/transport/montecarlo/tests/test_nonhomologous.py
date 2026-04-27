import pytest
from numpy.testing import assert_almost_equal

from tardis.transport.montecarlo.nonhomologous_grid import depressed_quartic


@pytest.mark.parametrize(
    ["A", "B", "C", "D", "E", "expected_roots"],
    [
        # x^4 - 14x^3 + 71x^2 - 154x + 120 = 0
        # roots 2, 3, 4, 5
        (1.0, -14.0, 71.0, -154.0, 120.0, [5.0, 4.0, 3.0, 2.0]),
        # x^4 - x^3 = 0
        # x^3(x + 1) = 0; only real root other than 0 is 1.0
        (1.0, -1.0, 0.0, 0.0, 0.0, [1.0, 0.0, 0.0, 0.0]),
        # x^4 - 10x^3 + 35x^2 - 50x + 24 = 0
        # roots 1, 2, 3, 4
        (1.0, -10.0, 35.0, -50.0, 24.0, [4.0, 3.0, 2.0, 1.0]),
    ],
)
def test_depressed_quartic(A, B, C, D, E, expected_roots):
    roots = depressed_quartic(A, B, C, D, E)
    assert_almost_equal(roots, expected_roots)
