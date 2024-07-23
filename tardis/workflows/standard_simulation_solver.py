import logging
from astropy import units as u

from tardis.simulation.convergence import ConvergenceSolver

# logging support
logger = logging.getLogger(__name__)


class StandardSimulationSolver:
    def __init__(self, convergence_strategy):
        self.simulation_state = None
        self.spectrum_solver = None
        self.transport_solver = None
        self.luminosity_requested = 0 * u.erg / u.s

        # Convergence
        self.convergence_strategy = convergence_strategy
        self.converged = False
        self.consecutive_converges_count = 0

        # Convergence solvers
        self.t_rad_convergence_solver = ConvergenceSolver(
            self.convergence_strategy.t_rad
        )
        self.w_convergence_solver = ConvergenceSolver(
            self.convergence_strategy.w
        )
        self.t_inner_convergence_solver = ConvergenceSolver(
            self.convergence_strategy.t_inner
        )

    def _get_convergence_status(
        self, t_rad, w, t_inner, estimated_t_rad, estimated_w, estimated_t_inner
    ):
        t_rad_converged = self.t_rad_convergence_solver.get_convergence_status(
            t_rad.value,
            estimated_t_rad.value,
            self.simulation_state.no_of_shells,
        )

        w_converged = self.w_convergence_solver.get_convergence_status(
            w, estimated_w, self.simulation_state.no_of_shells
        )

        t_inner_converged = (
            self.t_inner_convergence_solver.get_convergence_status(
                t_inner.value,
                estimated_t_inner.value,
                1,
            )
        )

        if np.all([t_rad_converged, w_converged, t_inner_converged]):
            hold_iterations = self.convergence_strategy.hold_iterations
            self.consecutive_converges_count += 1
            logger.info(
                f"Iteration converged {self.consecutive_converges_count:d}/{(hold_iterations + 1):d} consecutive "
                f"times."
            )
            # If an iteration has converged, require hold_iterations more
            # iterations to converge before we conclude that the Simulation
            # is converged.
            return self.consecutive_converges_count == hold_iterations + 1

        self.consecutive_converges_count = 0
        return False

    def check_convergence(self, emitted_luminosity):
        (
            estimated_t_rad,
            estimated_dilution_factor,
        ) = self.transport_solver.transport_state.calculate_radiationfield_properties()

        estimated_t_inner = self.estimate_t_inner(
            self.simulation_state.t_inner,
            self.luminosity_requested,
            emitted_luminosity,
            t_inner_update_exponent=self.convergence_strategy.t_inner_update_exponent,
        )

        converged = self._get_convergence_status(
            self.simulation_state.t_radiative,
            self.simulation_state.dilution_factor,
            self.simulation_state.t_inner,
            estimated_t_rad,
            estimated_dilution_factor,
            estimated_t_inner,
        )

        return converged

    def solve(self):
        while not converged:
            solve_plasma()
            run_montecarlo()
            converged = check_convergence()

        run_montecarlo(final=True)
