from tardis.plasma.equilibrium.rates.photoionization_strengths import (
    AnalyticCorrectedPhotoionizationCoeffSolver,
    EstimatedPhotoionizationCoeffSolver,
    SpontaneousRecombinationCoeffSolver,
)
from tardis.plasma.equilibrium.rates.util import (
    reindex_ionization_rate_dataframe,
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

    def solve(
        self,
        radiation_field,
        electron_energy_distribution,
        lte_level_population,
        level_population,
        lte_ion_population,
        ion_population,
        partition_function,
        level_boltzmann_factor,
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

        spontaneous_recombination_rate_coeff = (
            self.spontaneous_recombination_rate_coeff_solver.solve(
                electron_energy_distribution.temperature
            )
        )

        # TODO: Update for non-Hydrogenic species
        fractional_level_population = (
            level_boltzmann_factor / partition_function
        )

        # Lucy 2003 Eq 14
        level_to_ion_population_factor = lte_level_population.values / (
            lte_ion_population.values
            * electron_energy_distribution.number_density
        )

        # used to scale the photoionization rate because we keep the level population
        # fixed while we calculated the ion number density
        photoionization_rate = (
            photoionization_rate_coeff * fractional_level_population
        )

        # Lucy 2003 Eq 20
        spontaneous_recombination_rate = (
            spontaneous_recombination_rate_coeff
            * level_to_ion_population_factor
            * electron_energy_distribution.number_density
        )

        photoionization_rate = reindex_ionization_rate_dataframe(
            photoionization_rate, recombination=False
        )

        spontaneous_recombination_rate = reindex_ionization_rate_dataframe(
            spontaneous_recombination_rate, recombination=True
        )

        return photoionization_rate, spontaneous_recombination_rate


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
        estimators_continuum,
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
        estimators_continuum : EstimatorsContinuum
            Estimators of the continuum radiation field properties.
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
            estimators_continuum,
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

        photoionization_rate = reindex_ionization_rate_dataframe(
            photoionization_rate, recombination=False
        )

        recombination_rate = reindex_ionization_rate_dataframe(
            recombination_rate, recombination=True
        )

        return photoionization_rate, recombination_rate
