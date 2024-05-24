import numpy as np


class ConvergenceSolver:
    def __init__(self, strategy):
        """Convergence solver. Sets convergence strategy and assigns a method
        to the converge property.

        Parameters
        ----------
        strategy : string
            Convergence strategy for the physical property

        Raises
        ------
        NotImplementedError
            Custom convergence type specified
        ValueError
            Unknown convergence type specified
        """
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
        """Damped convergence solver

        Parameters
        ----------
        value : np.float64
            The current value of the physical property
        estimated_value : np.float64
            The estimated value of the physical property

        Returns
        -------
        np.float64
            The converged value
        """
        return value + self.damping_factor * (estimated_value - value)

    def get_convergence_status(self, value, estimated_value, no_of_cells):
        """Get the status of convergence for the physical property

        Parameters
        ----------
        value : np.float64, Quantity
            The current value of the physical property
        estimated_value : np.float64, Quantity
            The estimated value of the physical property
        no_of_cells : np.int64
            The number of cells to measure convergence over

        Returns
        -------
        bool
            True if convergence is reached
        """
        convergence = abs(value - estimated_value) / estimated_value

        fraction_converged = (
            np.count_nonzero(convergence < self.threshold) / no_of_cells
        )
        return fraction_converged > self.threshold
