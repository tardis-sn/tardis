import numpy as np
from astropy import units as u
from astropy.units import Quantity

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
        self.memory = self.convergence_strategy.get('memory', 5)

        #Initializing the history for Anderson method
        self.x_history = []
        self.g_history = []

        if self.convergence_strategy.type in ("damped"):
            self.converge = self.damped_converge
        elif self.convergence_strategy.type in ("anderson"):
            self.converge = self.anderson_converge
        elif self.convergence_strategy.type in ("custom"):
            raise NotImplementedError(
                "Convergence strategy type is custom; "
                "you need to implement your specific treatment!"
            )
        else:
            raise ValueError(
                f"Convergence strategy type is "
                f"not damped, anderson or custom "
                f"- input is {self.convergence_strategy.type}"
            )

    def anderson_converge(self, value, estimated_value):
        """
        Anderson acceleration method for convergence.

        Parameters
        ----------
        value : Quantity or float
            The current value of the physical property.
        estimated_value : Quantity or float
            The estimated value of the physical property.

        Returns
        -------
        Quantity
            The accelerated converged value with the same units as the input value.

        Raises
        ------
        LinAlgError
            If the QR decomposition or the linear system solution fails, the method 
            should fall back to returning the current estimated value without acceleration.
        """
        if not isinstance(value, Quantity):
            value = value * u.K  
        if not isinstance(estimated_value, Quantity):
            estimated_value = estimated_value * u.K

        original_unit = value.unit
        if value.unit != estimated_value.unit:
            value = value.to(estimated_value.unit)

        value = value.value
        estimated_value = estimated_value.value

        x_new = estimated_value
        g_new = x_new - value

        self.x_history.append(value) 
        self.g_history.append(g_new)

        if len(self.x_history) > 1:
            m = min(self.memory, len(self.x_history) - 1)
            G_k = np.column_stack([self.g_history[i] - self.g_history[i - 1] for i in range(-m, 0)])
            X_k = np.column_stack([self.x_history[i] - self.x_history[i - 1] for i in range(-m, 0)])

            g_new = np.array(g_new).reshape(-1, 1)
            try:
                Q, R = np.linalg.qr(G_k)
                gamma_k = np.linalg.solve(R, Q.T @ g_new)
                x_accel = self.x_history[-1] - (X_k @ gamma_k + G_k @ gamma_k).flatten()
                self.x_history[-1] = x_accel
                self.g_history[-1] = estimated_value - x_accel
                return (x_accel * original_unit)
            except np.linalg.LinAlgError:
                pass

        return x_new * original_unit 

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
