import numpy as np
import numpy.testing as npt
import pytest

from tardis.simulation.convergence import ConvergenceSolver


@pytest.mark.parametrize("damping_factor, threshold", [(0.5, 0.05), (0, 0)])
def test_convergence_solver_init_damped(
    t_rad_strategy, damping_factor, threshold
):
    t_rad_strategy.damping_factor = damping_factor
    t_rad_strategy.threshold = threshold
    solver = ConvergenceSolver(t_rad_strategy)
    assert solver.damping_factor == damping_factor
    assert solver.threshold == threshold
    assert solver.converge == solver.damped_converge


@pytest.mark.parametrize(
    "strategy_type, expected_exception",
    [("custom", NotImplementedError), ("invalid", ValueError)],
)
def test_convergence_solver_init_invalid(
    t_rad_strategy, strategy_type, expected_exception
):
    t_rad_strategy.type = strategy_type
    with pytest.raises(expected_exception):
        ConvergenceSolver(t_rad_strategy)


@pytest.mark.parametrize(
    "value, estimated_value, expected_converged_value",
    [(10.0, 20.0, 15.0), (5.0, 15.0, 10.0), (0.0, 10.0, 5.0)],
)
def test_damped_converge(
    t_rad_strategy, value, estimated_value, expected_converged_value
):
    solver = ConvergenceSolver(t_rad_strategy)
    converged_value = solver.damped_converge(
        np.float64(value), np.float64(estimated_value)
    )
    npt.assert_almost_equal(converged_value, expected_converged_value)


def test_get_convergence_status(t_rad_strategy):
    solver = ConvergenceSolver(t_rad_strategy)
    value = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    estimated_value = np.array([1.01, 2.02, 3.03], dtype=np.float64)
    no_of_cells = np.int64(3)
    is_converged = solver.get_convergence_status(
        value, estimated_value, no_of_cells
    )
    assert is_converged

    value = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    estimated_value = np.array([2.0, 3.0, 4.0], dtype=np.float64)
    is_converged = solver.get_convergence_status(
        value, estimated_value, no_of_cells
    )
    assert not is_converged
