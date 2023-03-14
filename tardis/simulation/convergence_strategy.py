import logging
import numpy as np
from IPython.display import display

# Adding logging support
logger = logging.getLogger(__name__)


class ConvergenceStrategy(object):
    def __init__(self, config, model):
        """
        Create a new ConvergenceStrategy instance from a Configuration object and Model object.

        Parameters
        ----------
        config : tardis.io.config_reader.Configuration
        model: tardis.model.Radial1DModel
        """
        # Open up the convergence_strategy dict in the config
        for key, value in config.montecarlo.convergence_strategy.items():
            setattr(self, key, value)
        self.model = model
        self._converged = (
            False  # Whether the convergence has truly converged or not
        )
        self._consecutive_converges_count = (
            0  # Counts the number of consecutive converges
        )

    @property
    def converged(self):
        return self._converged

    @property
    def should_stop(self):
        return self.stop_if_converged and self._converged

    def convergence_step(
        self,
        estimated_t_rad,
        estimated_w,
        estimated_t_inner,
        iterations_executed,
    ):
        """
        Performs a single convergence step

        Parameters
        ----------
        estimated_t_rad : astropy.units.Quantity
            The estimated value of t_rad
        estimated_w : np.ndarray
            The estimated value of w
        estimated_t_inner : astropy.units.Quantity
            The estimated value of t_inner
        iterations_executed : int
            Number of iterations elapsed in the simulation (to decide if t_inner should be updated)
        """
        # Compute convergence status
        self._converged = self._get_convergence_status(
            estimated_t_rad, estimated_w, estimated_t_inner
        )

        # Get new values
        if self.type in ("damped"):
            next_t_rad = self.damped_converge(
                self.model.t_rad,
                estimated_t_rad,
                self.t_rad.damping_constant,
            )
            next_w = self.damped_converge(
                self.model.w, estimated_w, self.w.damping_constant
            )
            if (iterations_executed + 1) % self.lock_t_inner_cycles == 0:
                next_t_inner = self.damped_converge(
                    self.model.t_inner,
                    estimated_t_inner,
                    self.t_inner.damping_constant,
                )
            else:
                next_t_inner = self.model.t_inner
        elif self.type in ("custom"):
            raise NotImplementedError(
                "Convergence strategy type is custom; "
                "you need to implement your specific treatment!"
            )
        else:
            raise ValueError(
                f"Convergence strategy type is "
                f"not damped or custom "
                f"- input is {self.type}"
            )

        return next_t_rad, next_w, next_t_inner

    @staticmethod
    def damped_converge(value, estimated_value, damping_factor):
        return value + damping_factor * (estimated_value - value)

    def _get_convergence_status(
        self, estimated_t_rad, estimated_w, estimated_t_inner
    ):
        """
        Returns whether the convergence has truly converged or not

        Parameters
        ----------
        estimated_t_rad : astropy.units.Quantity
            The estimated value of t_rad
        estimated_w : np.ndarray
            The estimated value of w
        estimated_t_inner : astropy.units.Quantity
            The estimated value of t_inner
        """
        no_of_shells = self.model.no_of_shells

        t_rad, w, t_inner = self.model.t_rad, self.model.w, self.model.t_inner

        convergence_t_rad = (
            abs(t_rad - estimated_t_rad) / estimated_t_rad
        ).value
        convergence_w = abs(w - estimated_w) / estimated_w
        convergence_t_inner = (
            abs(t_inner - estimated_t_inner) / estimated_t_inner
        ).value

        fraction_t_rad_converged = (
            np.count_nonzero(convergence_t_rad < self.t_rad.threshold)
            / no_of_shells
        )

        t_rad_converged = fraction_t_rad_converged > self.fraction

        fraction_w_converged = (
            np.count_nonzero(convergence_w < self.w.threshold) / no_of_shells
        )

        w_converged = fraction_w_converged > self.fraction

        t_inner_converged = convergence_t_inner < self.t_inner.threshold

        if np.all([t_rad_converged, w_converged, t_inner_converged]):
            hold_iterations = self.hold_iterations
            self._consecutive_converges_count += 1
            logger.info(
                f"Iteration converged {self._consecutive_converges_count:d}/{(hold_iterations + 1):d} consecutive "
                f"times."
            )
            # If an iteration has converged, require hold_iterations more
            # iterations to converge before we conclude that the Simulation
            # is converged.
            return self._consecutive_converges_count == hold_iterations + 1
        else:
            self._consecutive_converges_count = 0
            return False
