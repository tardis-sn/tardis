import numpy as np


class ConvergenceSolver:
    def __init__(self, strategy):
        self.convergence_strategy = strategy
        self.damping_factor = self.convergence_strategy.damping_constant
        self.threshold = self.convergence_strategy.threshold

        if self.convergence_strategy.type in ("damped"):
            self.converge = self.damped_converge
        elif self.convergence_strategy.type in ("custom"):
            raise NotImplementedError(
                "Convergence strategy type is custom; "
                "you need to implement your specific treatment!"
            )
        else:
            raise ValueError(
                f"Convergence strategy type is "
                f"not damped or custom "
                f"- input is {self.convergence_strategy.type}"
            )

    def damped_converge(self, value, estimated_value):
        return value + self.damping_factor * (estimated_value - value)

    def get_convergence_status(self, value, estimated_value, no_of_cells):
        convergence = abs(value - estimated_value) / estimated_value

        fraction_converged = (
            np.count_nonzero(convergence < self.threshold) / no_of_cells
        )
        return fraction_converged > self.threshold
