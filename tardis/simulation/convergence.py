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
        elif self.convergence_strategy.type in ("adaptive_damped"):
            self.lambda_min = 0.1
            self.lambda_max = 1.0
            self.lambda_step = 0.05
            self.residual_old = None 
            self.damping_factor = 0.5  
            self.converge = self.adaptive_damped
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
    
    def adaptive_damped(self, value, estimated_value):
        """Adaptive damped convergence solver

        Dynamically updates the damping factor by locally searching around the current value and selecting the step that minimizes the residual

        Parameters
        ----------
        value : np.float64
            The current value of the physical property
        estimated_value : np.float64
           The estimated value of the physical property

        Returns
        -------
        np.float64
            Updated value after applying adaptive damping

        Notes
        -----
        The damping factor is updated in place and constrained in the interval [lambda_min, lambda_max]
        """
        delta = self.lambda_step 
        base = self.damping_factor 

        candidates = [base]
        if base - delta >= self.lambda_min:
            candidates.append(base - delta)
        if base + delta <= self.lambda_max:
            candidates.append(base + delta)

        best_lambda = base
        best_residual = None

        for lam in candidates:
            x_new = value + lam * (estimated_value - value)
            res = np.mean(np.abs((estimated_value - x_new) / (estimated_value)))
            if best_residual is None or res < best_residual:
                best_residual = res
                best_lambda = lam

        self.damping_factor = best_lambda

        return value + best_lambda * (estimated_value - value)

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
