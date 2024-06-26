from pathlib import Path

import numpy as np
import numpy.testing as npt
import pytest

from tardis.io.configuration.config_reader import Configuration
from tardis.simulation.convergence import ConvergenceSolver


@pytest.fixture(scope="function")
def config(example_configuration_dir: Path):
    return Configuration.from_yaml(
        example_configuration_dir / "tardis_configv1_verysimple.yml"
    )


@pytest.fixture(scope="function")
def strategy(config):
    return config.montecarlo.convergence_strategy.t_rad


def test_convergence_solver_init_damped(strategy):
    solver = ConvergenceSolver(strategy)
    assert solver.damping_factor == 0.5
    assert solver.threshold == 0.05
    assert solver.converge == solver.damped_converge


def test_convergence_solver_init_custom(strategy):
    strategy.type = "custom"
    with pytest.raises(NotImplementedError):
        ConvergenceSolver(strategy)


def test_convergence_solver_init_invalid(strategy):
    strategy.type = "invalid"
    with pytest.raises(ValueError):
        ConvergenceSolver(strategy)


def test_damped_converge(strategy):
    solver = ConvergenceSolver(strategy)
    value = np.float64(10.0)
    estimated_value = np.float64(20.0)
    converged_value = solver.damped_converge(value, estimated_value)
    npt.assert_almost_equal(converged_value, 15.0)


def test_get_convergence_status(strategy):
    solver = ConvergenceSolver(strategy)
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
