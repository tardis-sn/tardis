import numpy as np
import pandas as pd

from tardis import constants as const
from tardis.transport.montecarlo.estimators.util import (
    bound_free_estimator_array2frame,
    integrate_array_by_blocks,
)

C = const.c.cgs.value
H = const.h.cgs.value
K_B = const.k_B.cgs.value


class SpontaneousRecombinationCoeffSolver:
    def __init__(
        self,
        photoionization_cross_sections,
    ):
        self.photoionization_cross_sections = photoionization_cross_sections
        self.nu = self.photoionization_cross_sections.nu.values

    @property
    def common_prefactor(self):
        return (
            4.0
            * np.pi
            * self.photoionization_cross_sections.x_sect
            / (H * self.nu)
        )

    def calculate_photoionization_boltzmann_factor(self, electron_temperature):
        return np.exp(-self.nu[np.newaxis].T / electron_temperature * (H / K_B))

    def solve(self, electron_temperature):
        """
        Calculate the spontaneous recombination rate coefficient.

        Parameters
        ----------
        electron_temperature : u.Quantity
            Electron temperature in each cell.

        Returns
        -------
        pd.DataFrame
            The calculated spontaneous recombination rate coefficient.

        Notes
        -----
        Equation 13 in Lucy 2003.
        """
        # need to fix array multiplication
        prefactor = self.common_prefactor * (2 * H * self.nu**3.0) / (C**2.0)
        photoionization_boltzmann_factor = pd.DataFrame(
            self.calculate_photoionization_boltzmann_factor(
                electron_temperature
            ),
            index=prefactor.index,
        )
        spontaneous_recombination_rate_coeff = (
            photoionization_boltzmann_factor.multiply(
                prefactor,
                axis=0,
            )
        )
        return spontaneous_recombination_rate_coeff


class AnalyticPhotoionizationCoeffSolver(SpontaneousRecombinationCoeffSolver):
    def __init__(
        self,
        photoionization_cross_sections,
    ):
        super().__init__(photoionization_cross_sections)

        self.photoionization_block_references = np.pad(
            self.photoionization_cross_sections.nu.groupby(level=[0, 1, 2])
            .count()
            .values.cumsum(),
            [1, 0],
        )

        self.photoionization_index = (
            self.photoionization_cross_sections.index.unique()
        )

    def calculate_mean_intensity_photoionization_df(
        self,
        dilute_blackbody_radiationfield_state,
    ):
        mean_intensity = (
            dilute_blackbody_radiationfield_state.calculate_mean_intensity(
                self.nu
            )
        )
        return pd.DataFrame(
            mean_intensity,
            index=self.photoionization_cross_sections.index,
            columns=np.arange(
                len(dilute_blackbody_radiationfield_state.temperature)
            ),
        )

    def calculate_photoionization_rate_coeff(
        self,
        mean_intensity_photoionization_df,
    ):
        """
        Calculate the photoionization rate coefficient.

        Parameters
        ----------
        dilute_blackbody_radiationfield_state : DiluteBlackBodyRadiationFieldState
            A dilute black body radiation field state.

        Returns
        -------
        pd.DataFrame
            The calculated photoionization rate coefficient.

        Notes
        -----
        Equation 16 in Lucy 2003.
        """
        photoionization_rate_coeff = mean_intensity_photoionization_df.multiply(
            self.common_prefactor,
            axis=0,
        )
        photoionization_rate_coeff = integrate_array_by_blocks(
            photoionization_rate_coeff.values,
            self.nu,
            self.photoionization_block_references,
        )
        photoionization_rate_coeff = pd.DataFrame(
            photoionization_rate_coeff,
            index=self.photoionization_index,
        )
        return photoionization_rate_coeff

    def calculate_stimulated_recombination_rate_coeff(
        self,
        mean_intensity_photoionization_df,
        photoionization_boltzmann_factor,
    ):
        """
        Calculate the photoionization rate coefficient.

        Parameters
        ----------
        mean_intensity_photoionization_df : pd.DataFrame
            Mean intensity at each photoionization frequency.
        photoionization_boltzmann_factor : np.ndarray
            Boltzmann factor for each photoionization frequency.

        Returns
        -------
        pd.DataFrame
            The stimulated recombination rate coefficient.

        Notes
        -----
        Equation 15 in Lucy 2003.
        """
        stimulated_recombination_rate_coeff = (
            mean_intensity_photoionization_df * photoionization_boltzmann_factor
        )

        stimulated_recombination_rate_coeff = (
            mean_intensity_photoionization_df.multiply(
                self.common_prefactor,
                axis=0,
            )
        )
        stimulated_recombination_rate_coeff = integrate_array_by_blocks(
            stimulated_recombination_rate_coeff.values,
            self.nu,
            self.photoionization_block_references,
        )
        stimulated_recombination_rate_coeff = pd.DataFrame(
            stimulated_recombination_rate_coeff,
            index=self.photoionization_index,
        )
        return stimulated_recombination_rate_coeff

    def solve(
        self,
        dilute_blackbody_radiationfield_state,
        electron_temperature,
    ):
        """
        Prepares the ionization and recombination coefficients by grouping them for
        ion numbers.

        Parameters
        ----------
        dilute_blackbody_radiationfield_state : DiluteBlackBodyRadiationFieldState
            The dilute black body radiation field state.
        electron_temperature : u.Quantity
            Electron temperature in each shell.

        Returns
        -------
        photoionization_rate_coeff
            Photoionization rate coefficient grouped by atomic number and ion number.
        recombination_rate_coeff
            Radiative recombination rate coefficient grouped by atomic number and ion number.
        """
        photoionization_boltzmann_factor = (
            self.calculate_photoionization_boltzmann_factor(
                electron_temperature
            )
        )

        mean_intensity_photoionization_df = (
            self.calculate_mean_intensity_photoionization_df(
                dilute_blackbody_radiationfield_state
            )
        )

        photoionization_rate_coeff = self.calculate_photoionization_rate_coeff(
            mean_intensity_photoionization_df,
        )
        # Equation 15 Lucy 2003. Must be multiplied by Saha LTE factor Phi_ik
        stimulated_recombination_rate_coeff = (
            self.calculate_stimulated_recombination_rate_coeff(
                mean_intensity_photoionization_df,
                photoionization_boltzmann_factor,
            )
        )

        return (
            photoionization_rate_coeff,
            stimulated_recombination_rate_coeff,
        )


class EstimatedPhotoionizationCoeffSolver:
    def __init__(
        self,
        level2continuum_edge_idx,
    ):
        self.level2continuum_edge_idx = level2continuum_edge_idx

    def solve(
        self,
        radfield_mc_estimators,
        time_simulation,
        volume,
    ):
        """
        Solve for the continuum properties.

        Parameters
        ----------
        radfield_mc_estimators : RadiationFieldMCEstimators
            The Monte Carlo estimators for the radiation field.
        time_simulation : float
            The simulation time.
        volume : float
            The volume of the cells.

        Returns
        -------
        ContinuumProperties
            The calculated continuum properties.
        """
        photoionization_normalization = (time_simulation * volume * H) ** -1

        photoionization_rate_coeff = bound_free_estimator_array2frame(
            radfield_mc_estimators.photo_ion_estimator,
            self.level2continuum_edge_idx,
        )
        photoionization_rate_coeff *= photoionization_normalization

        stimulated_recombination_rate_coeff = bound_free_estimator_array2frame(
            radfield_mc_estimators.stim_recomb_estimator,
            self.level2continuum_edge_idx,
        )
        stimulated_recombination_rate_coeff *= photoionization_normalization

        return photoionization_rate_coeff, stimulated_recombination_rate_coeff
