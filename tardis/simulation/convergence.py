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
            self.damping_factor = 0.5  # starting point
            self.converge = self.adaptive_damped
        elif self.convergence_strategy.type in ("RobbinsMonro"):
            self.iteration = 1
            self.converge = self.harmonic_damped
        elif self.convergence_strategy.type in ("custom"):
            raise NotImplementedError(
                "Convergence strategy type is custom; "
                "you need to implement your specific treatment!"
            )
        else:
            raise ValueError(
                f"Convergence strategy type is "
                f"not damped, adaptive damped, Robbins-Monro, or custom "
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
    
    # def adaptive_damped(self, value, estimated_value):
    #     """
    #     Monotonic non‐decreasing adaptive damping:
    #     - Start at self.damping_factor 
    #     - Candidate lambdas: {λ, λ + Δ}, clipped to [λ, λ_max]
    #     - Compute residual r = mean(|estimated - x_new| / |estimated|)
    #     - Pick the candidate that minimizes r
    #     - Update self.damping_factor to that candidate (never decreases)
    #     - Return the damped update with that λ
    #     """
    #     base  = self.damping_factor
    #     delta = self.lambda_step

    #     # Only allowing lambda to stay the same or increase
    #     candidates = [base]
    #     up = min(base + delta, self.lambda_max)
    #     if up > base:
    #         candidates.append(up)

    #     best_lambda   = base
    #     best_residual = None

    #     for lam in candidates:
    #         x_new = value + lam * (estimated_value - value)
    #         # normalized residual
    #         res = np.mean(np.abs((estimated_value - x_new) / estimated_value))
    #         if best_residual is None or res < best_residual:
    #             best_residual = res
    #             best_lambda   = lam

    #     # This assignment can only increase or hold λ
    #     self.damping_factor = best_lambda

    #     return value + best_lambda * (estimated_value - value)

    # def adaptive_damped(self, value, estimated_value):
    #     """
    #     Monotonically non-increasing adaptive damping: λ₀ = 1.0
    #     At each iteration n:
    #       rₙ = mean(|est - val|/|est|)
    #       if rₙ >= rₙ₋₁: λₙ₊₁ = max(λₙ - Δ, λ_min)
    #       else:         λₙ₊₁ = λₙ
    #     """
    #     # Compute the new normalized residual
    #     resid_new = np.mean(np.abs((estimated_value - value) / estimated_value))

    #     if self.residual_old is not None:
    #         # r_n is worse or unchanged?
    #         if resid_new >= self.residual_old:
    #             # step down
    #             self.damping_factor = max(self.damping_factor - self.lambda_step,
    #                                       self.lambda_min)
    #         # else: leave damping_factor unchanged

    #     # store for next iteration
    #     self.residual_old = resid_new

    #     # apply update
    #     return value + self.damping_factor * (estimated_value - value)


    def adaptive_damped(self, value, estimated_value):
        """
        Split(Local-search) damping:
        - Start at 0.5
        - Form candidates {λ, λ-Δ, λ+Δ} clipped to [0,1]
        - Compute the residual for each: r = mean(|x* - (x + λ*(x*-x))|)
        - Pick the λ that minimizes r, update self.damping_factor
        - Return the damped update with that λ
        """
        delta = self.lambda_step # 0.05
        base = self.damping_factor # 0.5

        candidates = [base]
        if base - delta >= self.lambda_min:
            candidates.append(base - delta)
        if base + delta <= self.lambda_max:
            candidates.append(base + delta)

        best_lambda = base
        best_residual = None

        for lam in candidates:
            x_new = value + lam * (estimated_value - value)
            res = np.mean(abs(np.abs((estimated_value - x_new) / (estimated_value))))
            if best_residual is None or res < best_residual:
                best_residual = res
                best_lambda = lam

        self.damping_factor = best_lambda

        return value + best_lambda * (estimated_value - value)
    
    def harmonic_damped(self, value, estimated_value):
        """
        Robbins-Monro convergence algorithm: λ_n = 1 / n
        Damping decreases with each iteration.
        """
        self.damping_factor = 1.0 / self.iteration
        self.iteration += 1 

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
