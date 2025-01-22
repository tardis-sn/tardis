from tardis.plasma.equilibrium.rates.photoionization_strengths import (
    AnalyticPhotoionizationCoeffSolver,
    EstimatedPhotoionizationCoeffSolver,
    SpontaneousRecombinationCoeffSolver,
)


class AnalyticPhotoionizationRateSolver:
    """Solve the photoionization and spontaneous recombination rates in the
    case where the radiation field is computed analytically.
    """

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
        level_number_density,
        ion_number_density,
        electron_number_density,
        saha_factor,
    ):
        """Compute the photoionization and spontaneous recombination rates

        Parameters
        ----------
        photoionization_rate_coeff : pd.DataFrame
            The photoionization rate coefficients for each transition.
            Columns are cells.
        stimulated_recombination_rate_coeff : pd.DataFrame
            The stimulated recombination rate coefficients for each transition.
            Columns are cells.
        spontaneous_recombination_rate_coeff : pd.DataFrame
            The spontaneous recombination rate coefficients for each transition.
            Columns are cells.
        level_number_density : pd.DataFrame
            The electron energy level number density. Columns are cells.
        ion_number_density : pd.DataFrame
            The ion number density. Columns are cells.
        electron_number_density : u.Quantity
            The free electron number density per cell.
        saha_factor : pd.DataFrame
            The LTE population factor. Columns are cells.

        Returns
        -------
        pd.DataFrame
            Photoionization rate for each electron energy level. Columns are cells
        pd.DataFrame
            Spontaneous recombination rate for each electron energy level. Columns are cells
        """
        photoionization_rate = (
            photoionization_rate_coeff * level_number_density
            - saha_factor
            * stimulated_recombination_rate_coeff
            * ion_number_density
            * electron_number_density
        )
        spontaneous_recombination_rate = (
            saha_factor
            * spontaneous_recombination_rate_coeff
            * ion_number_density
            * electron_number_density
        )

        return photoionization_rate, spontaneous_recombination_rate

    def solve(
        self,
        radiation_field,
        electron_energy_distribution,
        level_number_density,
        ion_number_density,
        saha_factor,
    ):
        """Solve the photoionization and spontaneous recombination rates in the
        case where the radiation field is not estimated.

        Parameters
        ----------
        radiation_field : RadiationField
            A radiation field that can compute its mean intensity.
        electron_energy_distribution : ThermalElectronEnergyDistribution
            Electron properties.
        level_number_density : pd.DataFrame
            Electron energy level number density. Columns are cells.
        ion_number_density : pd.DataFrame
            Ion number density. Columns are cells.
        saha_factor : pd.DataFrame
            Saha factor: the LTE level number density divided by the LTE ion
            number density and the electron number density.

        Returns
        -------
        pd.DataFrame
            Photoionization rate. Columns are cells.
        pd.DataFrame
            Spontaneous recombination rate. Columns are cells.
        """
        photoionization_rate_coeff_solver = AnalyticPhotoionizationCoeffSolver(
            self.photoionization_cross_sections
        )

        photoionization_rate_coeff, stimulated_recombination_rate_coeff = (
            photoionization_rate_coeff_solver.solve(
                radiation_field,
                electron_energy_distribution.temperature,
            )
        )

        spontaneous_recombination_rate_coeff = (
            self.spontaneous_recombination_rate_coeff_solver.solve(
                electron_energy_distribution.temperature
            )
        )

        return self.compute_rates(
            photoionization_rate_coeff,
            stimulated_recombination_rate_coeff,
            spontaneous_recombination_rate_coeff,
            level_number_density,
            ion_number_density,
            electron_energy_distribution.number_density,
            saha_factor,
        )


class EstimatedPhotoionizationRateSolver(AnalyticPhotoionizationRateSolver):
    """Solve the photoionization and spontaneous recombination rates in the
    case where the radiation field is estimated by Monte Carlo processes.
    """

    def __init__(
        self, photoionization_cross_sections, level2continuum_edge_idx
    ):
        super().__init__(
            photoionization_cross_sections,
        )
        self.level2continuum_edge_idx = level2continuum_edge_idx

    def solve(
        self,
        electron_energy_distribution,
        radfield_mc_estimators,
        time_simulation,
        volume,
        level_number_density,
        ion_number_density,
        saha_factor,
    ):
        """Solve the photoionization and spontaneous recombination rates in the
        case where the radiation field is estimated by Monte Carlo processes.

        Parameters
        ----------
        electron_energy_distribution : ThermalElectronEnergyDistribution
            Electron properties.
        radfield_mc_estimators : RadiationFieldMCEstimators
            Estimators of the radiation field properties.
        time_simulation : u.Quantity
            Time of simulation.
        volume : u.Quantity
            Volume per cell.
        level_number_density : pd.DataFrame
            Electron energy level number density. Columns are cells.
        ion_number_density : pd.DataFrame
            Ion number density. Columns are cells.
        saha_factor : pd.DataFrame
            Saha factor: the LTE level number density divided by the LTE ion
            number density and the electron number density.

        Returns
        -------
        pd.DataFrame
            Photoionization rate. Columns are cells.
        pd.DataFrame
            Spontaneous recombination rate. Columns are cells.
        """
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
                electron_energy_distribution.temperature
            )
        )

        return self.compute_rates(
            photoionization_rate_coeff,
            stimulated_recombination_rate_coeff,
            spontaneous_recombination_rate_coeff,
            level_number_density,
            ion_number_density,
            electron_energy_distribution.number_density,
            saha_factor,
        )
