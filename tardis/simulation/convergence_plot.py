import logging
import numpy as np

logger = logging.getLogger(__name__)


class ConvergencePlot:
    @staticmethod
    def damped_converge(self,value, estimated_value, damping_factor):
        return value + damping_factor * (estimated_value - value)

    # @staticmethod
    # def _get_convergence_status(
    #     t_rad, w, t_inner, estimated_t_rad, estimated_w, estimated_t_inner,
    #     no_of_shells, convergence_strategy
    # ):
    #     convergence_t_rad = (
    #         abs(t_rad - estimated_t_rad) / estimated_t_rad
    #     ).value
    #     convergence_w = abs(w - estimated_w) / estimated_w
    #     convergence_t_inner = (
    #         abs(t_inner - estimated_t_inner) / estimated_t_inner
    #     ).value

    #     fraction_t_rad_converged = (
    #         np.count_nonzero(
    #             convergence_t_rad < convergence_strategy.t_rad.threshold
    #         )
    #         / no_of_shells
    #     )

    #     t_rad_converged = (
    #         fraction_t_rad_converged > convergence_strategy.fraction
    #     )

    #     fraction_w_converged = (
    #         np.count_nonzero(
    #             convergence_w < convergence_strategy.w.threshold
    #         )
    #         / no_of_shells
    #     )

    #     w_converged = fraction_w_converged > convergence_strategy.fraction

    #     t_inner_converged = (
    #         convergence_t_inner < convergence_strategy.t_inner.threshold
    #     )

    #     if np.all([t_rad_converged, w_converged, t_inner_converged]):
    #         return True
    #     else:
    #         return False

    def _get_convergence_status(
        self, t_rad, w, t_inner, estimated_t_rad, estimated_w, estimated_t_inner,
    ):
        no_of_shells = self.model.no_of_shells

        convergence_t_rad = (
            abs(t_rad - estimated_t_rad) / estimated_t_rad
        ).value
        convergence_w = abs(w - estimated_w) / estimated_w
        convergence_t_inner = (
            abs(t_inner - estimated_t_inner) / estimated_t_inner
        ).value

        fraction_t_rad_converged = (
            np.count_nonzero(
                convergence_t_rad < self.convergence_strategy.t_rad.threshold
            )
            / no_of_shells
        )

        t_rad_converged = (
            fraction_t_rad_converged > self.convergence_strategy.fraction
        )

        fraction_w_converged = (
            np.count_nonzero(
                convergence_w < self.convergence_strategy.w.threshold
            )
            / no_of_shells
        )

        w_converged = fraction_w_converged > self.convergence_strategy.fraction

        t_inner_converged = (
            convergence_t_inner < self.convergence_strategy.t_inner.threshold
        )

        if np.all([t_rad_converged, w_converged, t_inner_converged]):
            hold_iterations = self.convergence_strategy.hold_iterations
            self.consecutive_converges_count += 1
            logger.info(
                f"Iteration converged {self.consecutive_converges_count:d}/{(hold_iterations + 1):d} consecutive "
                f"times."
            )
            return self.consecutive_converges_count == hold_iterations + 1
        else:
            self.consecutive_converges_count = 0
            return False
