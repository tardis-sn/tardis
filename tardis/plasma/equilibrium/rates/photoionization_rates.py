from tardis.plasma.equilibrium.rates.photoionization_strengths import (
    AnalyticPhotoionizationCoeffSolver,
    EstimatedPhotoionizationCoeffSolver,
    SpontaneousRecombinationCoeffSolver,
)


class PhotoionizationRateSolver:
    def __init__(
        self, photoionization_cross_sections, level2continuum_edge_idx
    ):
        self.photoionization_cross_sections = photoionization_cross_sections
        self.level2continuum_edge_idx = level2continuum_edge_idx

    def solve(
        self, electron_temperature, n_i, n_k, n_e, saha_factor, type="analytic"
    ):
        if type == "analytic":
            # TODO: try merging the analytic and estimator based approaches
            # look at the Lucy 2003 equations and our MC estimator classes
            photoionization_rate_coeff_solver = (
                AnalyticPhotoionizationCoeffSolver(
                    self.photoionization_cross_sections
                )
            )
        elif type == "estimated":
            photoionization_rate_coeff_solver = (
                EstimatedPhotoionizationCoeffSolver(
                    self.level2continuum_edge_idx
                )
            )
        else:
            raise ValueError(f"Type {type} not supported")

        spontaneous_recombination_rate_coeff_solver = (
            SpontaneousRecombinationCoeffSolver(
                self.photoionization_cross_sections
            )
        )

        photoionization_rate_coeff, stimulated_recombination_rate_coeff = (
            # TODO: bifurcation of classes here is a problem
            photoionization_rate_coeff_solver.solve()
        )

        spontaneous_recombination_rate_coeff = (
            spontaneous_recombination_rate_coeff_solver.solve(
                electron_temperature
            )
        )

        # TODO: decide if these n_i, n_k, n_e numbers should be here (probably not)
        photoionization_rate = (
            photoionization_rate_coeff * n_i
            - saha_factor * stimulated_recombination_rate_coeff * n_k * n_e
        )
        spontaneous_recombination_rate = (
            saha_factor * spontaneous_recombination_rate_coeff * n_k * n_e
        )

        return photoionization_rate, spontaneous_recombination_rate
