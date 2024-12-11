from tardis.plasma.equilibrium.rates.photoionization_strengths import (
    AnalyticPhotoionizationCoeffSolver,
    EstimatedPhotoionizationCoeffSolver,
    SpontaneousRecombinationCoeffSolver,
)


class PhotoionizationRateSolver:
    def __init__():
        pass

    def solve(n_i, n_k, n_e, saha_factor, type="analytic"):
        if type == "analytic":
            photoionization_rate_coeff_solver = (
                AnalyticPhotoionizationCoeffSolver()
            )
        elif type == "estimated":
            photoionization_rate_coeff_solver = (
                EstimatedPhotoionizationCoeffSolver()
            )
        else:
            raise ValueError(f"Type {type} not supported")

        spontaneous_recombination_rate_coeff_solver = (
            SpontaneousRecombinationCoeffSolver()
        )

        photoionization_rate_coeff, stimulated_recombination_rate_coeff = (
            photoionization_rate_coeff_solver.solve()
        )

        spontaneous_recombination_rate_coeff = (
            spontaneous_recombination_rate_coeff_solver.solve()
        )

        photoionization_rate = (
            photoionization_rate_coeff * n_i
            - saha_factor * stimulated_recombination_rate_coeff * n_k * n_e
        )
        spontaneous_recombination_rate = (
            saha_factor * spontaneous_recombination_rate_coeff * n_k * n_e
        )

        return photoionization_rate, spontaneous_recombination_rate
