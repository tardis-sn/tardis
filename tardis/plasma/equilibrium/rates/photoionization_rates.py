from tardis.plasma.equilibrium.rates.photoionization_strengths import (
    AnalyticPhotoionizationCoeffSolver,
    EstimatedPhotoionizationCoeffSolver,
    SpontaneousRecombinationCoeffSolver,
)


class AnalyticPhotoionizationRateSolver:
    def __init__(self, photoionization_cross_sections):
        self.photoionization_cross_sections = photoionization_cross_sections

        self.spontaneous_recombination_rate_coeff_solver = (
            SpontaneousRecombinationCoeffSolver(
                self.photoionization_cross_sections
            )
        )

    def compute_rates(
        self,
        photoionization_rate_coeff,
        stimulated_recombination_rate_coeff,
        spontaneous_recombination_rate_coeff,
        n_i,
        n_k,
        n_e,
        saha_factor,
    ):
        # TODO: decide if these n_i, n_k, n_e numbers should be here (probably not)
        photoionization_rate = (
            photoionization_rate_coeff * n_i
            - saha_factor * stimulated_recombination_rate_coeff * n_k * n_e
        )
        spontaneous_recombination_rate = (
            saha_factor * spontaneous_recombination_rate_coeff * n_k * n_e
        )

        return photoionization_rate, spontaneous_recombination_rate

    def solve(self, electron_temperature, n_i, n_k, n_e, saha_factor):
        photoionization_rate_coeff_solver = AnalyticPhotoionizationCoeffSolver(
            self.photoionization_cross_sections
        )

        photoionization_rate_coeff, stimulated_recombination_rate_coeff = (
            photoionization_rate_coeff_solver.solve()
        )

        spontaneous_recombination_rate_coeff = (
            self.spontaneous_recombination_rate_coeff_solver.solve(
                electron_temperature
            )
        )

        return self.compute_rates(
            photoionization_rate_coeff,
            stimulated_recombination_rate_coeff,
            spontaneous_recombination_rate_coeff,
            n_i,
            n_k,
            n_e,
            saha_factor,
        )


class EstimatedPhotoionizationRateSolver(AnalyticPhotoionizationRateSolver):
    def __init__(
        self, photoionization_cross_sections, level2continuum_edge_idx
    ):
        super().__init__(
            photoionization_cross_sections,
        )
        self.level2continuum_edge_idx = level2continuum_edge_idx

    def solve(
        self,
        electron_temperature,
        radfield_mc_estimators,
        time_simulation,
        volume,
        n_i,
        n_k,
        n_e,
        saha_factor,
    ):
        photoionization_rate_coeff_solver = EstimatedPhotoionizationCoeffSolver(
            self.level2continuum_edge_idx
        )

        photoionization_rate_coeff, stimulated_recombination_rate_coeff = (
            photoionization_rate_coeff_solver.solve(
                radfield_mc_estimators,
                time_simulation,
                volume,
            )
        )

        spontaneous_recombination_rate_coeff = (
            self.spontaneous_recombination_rate_coeff_solver.solve(
                electron_temperature
            )
        )

        return self.compute_rates(
            photoionization_rate_coeff,
            stimulated_recombination_rate_coeff,
            spontaneous_recombination_rate_coeff,
            n_i,
            n_k,
            n_e,
            saha_factor,
        )
