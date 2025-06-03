from tardis.plasma.equilibrium.rates.photoionization_strengths import (
    AnalyticCorrectedPhotoionizationCoeffSolver,
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

    @staticmethod
    def __reindex_ionization_rate_dataframe(
        rate_dataframe, recombination=False
    ):
        rate_dataframe.index.names = [
            "atomic_number",
            "ion_number",
            "level_number_source",
        ]

        rate_dataframe = rate_dataframe.reset_index()

        if recombination:
            rate_dataframe["ion_number_destination"] = rate_dataframe[
                "ion_number"
            ]
            rate_dataframe["ion_number_source"] = (
                rate_dataframe["ion_number"] + 1
            )
        else:
            rate_dataframe["ion_number_source"] = rate_dataframe["ion_number"]
            rate_dataframe["ion_number_destination"] = (
                rate_dataframe["ion_number"] + 1
            )

        # ionized electrons are assumed to leave the ion in the ground state for now
        rate_dataframe["level_number_destination"] = 0

        not_fully_ionized_mask = (
            rate_dataframe["atomic_number"] != rate_dataframe["ion_number"]
        )

        rate_dataframe = rate_dataframe[not_fully_ionized_mask]

        rate_dataframe = rate_dataframe.set_index(
            [
                "atomic_number",
                "ion_number",
                "ion_number_source",
                "ion_number_destination",
                "level_number_source",
                "level_number_destination",
            ]
        )

        return rate_dataframe

    def solve(
        self,
        radiation_field,
        electron_energy_distribution,
        lte_level_population,
        level_population,
        lte_ion_population,
        ion_population,
    ):
        """Solve the photoionization and spontaneous recombination rates in the
        case where the radiation field is not estimated.

        Parameters
        ----------
        radiation_field : RadiationField
            A radiation field that can compute its mean intensity.
        electron_energy_distribution : ThermalElectronEnergyDistribution
            Electron properties.
        lte_level_population : pd.DataFrame
            LTE level number density. Columns are cells.
        level_population : pd.DataFrame
            Estimated level number density. Columns are cells.
        lte_ion_population : pd.DataFrame
            LTE ion number density. Columns are cells.
        ion_population : pd.DataFrame
            Estimated ion number density. Columns are cells.

        Returns
        -------
        pd.DataFrame
            Photoionization rate. Columns are cells.
        pd.DataFrame
            Spontaneous recombination rate. Columns are cells.
        """
        photoionization_rate_coeff_solver = (
            AnalyticCorrectedPhotoionizationCoeffSolver(
                self.photoionization_cross_sections
            )
        )

        photoionization_rate_coeff = photoionization_rate_coeff_solver.solve(
            radiation_field,
            electron_energy_distribution.temperature,
            lte_level_population,
            level_population,
            lte_ion_population,
            ion_population,
        )

        recombination_rate_coeff = (
            self.spontaneous_recombination_rate_coeff_solver.solve(
                electron_energy_distribution.temperature
            )
        )

        photoionization_rate = photoionization_rate_coeff * level_population

        recombination_rate = (
            recombination_rate_coeff
            * level_population
            * electron_energy_distribution.number_density
        )

        photoionization_rate = self.__reindex_ionization_rate_dataframe(
            photoionization_rate, recombination=False
        )

        recombination_rate = self.__reindex_ionization_rate_dataframe(
            recombination_rate, recombination=True
        )

        return photoionization_rate, recombination_rate


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
        level_population,
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
        level_population : pd.DataFrame
            Electron energy level number density. Columns are cells.

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

        photoionization_rate_coeff = photoionization_rate_coeff_solver.solve(
            radfield_mc_estimators,
            time_simulation,
            volume,
        )

        spontaneous_recombination_rate_coeff = (
            self.spontaneous_recombination_rate_coeff_solver.solve(
                electron_energy_distribution.temperature
            )
        )

        photoionization_rate = photoionization_rate_coeff * level_population

        recombination_rate = (
            spontaneous_recombination_rate_coeff
            * level_population
            * electron_energy_distribution.number_density
        )

        photoionization_rate = self.__reindex_ionization_rate_dataframe(
            photoionization_rate, recombination=False
        )

        recombination_rate = self.__reindex_ionization_rate_dataframe(
            recombination_rate, recombination=True
        )

        return photoionization_rate, recombination_rate
