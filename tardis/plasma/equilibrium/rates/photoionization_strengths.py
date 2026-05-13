import astropy.units as u
import numpy as np
import pandas as pd

from tardis import constants as const
from tardis.transport.montecarlo.estimators.util import (
    bound_free_estimator_array2frame,
    integrate_array_by_blocks,
)

C = const.c.cgs
H = const.h.cgs
K_B = const.k_B.cgs


class SpontaneousRecombinationCoeffSolver:
    def __init__(
        self,
        photoionization_cross_sections,
    ):
        self.photoionization_cross_sections = photoionization_cross_sections
        self.nu = self.photoionization_cross_sections.nu.values * u.Hz

        self.photoionization_block_references = np.pad(
            self.photoionization_cross_sections.nu.groupby(level=[0, 1, 2])
            .count()
            .values.cumsum(),
            [1, 0],
        )

        self.photoionization_index = (
            self.photoionization_cross_sections.index.unique()
        )

    @property
    def common_prefactor(self):
        """Used to multiply with both spontaneous recombination and
        photoionization coefficients. Lucy 2003 Eq 13, 15, 16.

        Returns
        -------
        pd.DataFrame
            A dataframe of the prefactor.
        """
        return (
            4.0
            * np.pi
            * self.photoionization_cross_sections.x_sect
            / (H * self.nu)
        )

    def calculate_photoionization_boltzmann_factor(self, electron_temperature):
        """Calculate the Boltzmann factor at each photoionization frequency

        Parameters
        ----------
        electron_temperature : Quantity
            Electron temperature in each shell.

        Returns
        -------
        numpy.ndarray
            The Boltzmann factor per shell per photoionization frequency.
        """
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
        Equation 13 in Lucy 2003, missing the factor from Eq 14.
        """
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
        spontaneous_recombination_rate_coeff_integrated = (
            integrate_array_by_blocks(
                spontaneous_recombination_rate_coeff.to_numpy(),
                self.nu.value,
                self.photoionization_block_references,
            )
        )

        spontaneous_recombination_rate_coeff_df = pd.DataFrame(
            spontaneous_recombination_rate_coeff_integrated,
            index=self.photoionization_index,
        )

        # Lymann continuum handling
        spontaneous_recombination_rate_coeff_df.loc[(1, 0, 0)] = 0.0

        return spontaneous_recombination_rate_coeff_df


class AnalyticPhotoionizationCoeffSolver(SpontaneousRecombinationCoeffSolver):
    def __init__(
        self,
        photoionization_cross_sections,
    ):
        super().__init__(photoionization_cross_sections)

    def calculate_mean_intensity_photoionization_df(
        self,
        dilute_blackbody_radiationfield_state,
    ):
        """Calculates the mean intensity of the radiation field at each photoionization frequency.

        Parameters
        ----------
        dilute_blackbody_radiationfield_state : DilutePlanckianRadiationField
            The radiation field.

        Returns
        -------
        pd.DataFrame
            DataFrame of mean intensities indexed by photoionization levels and
            columns of cells.
        """
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
            stimulated_recombination_rate_coeff.multiply(
                self.common_prefactor,
                axis=0,
            )
        )
        stimulated_recombination_rate_coeff = integrate_array_by_blocks(
            stimulated_recombination_rate_coeff.values,
            self.nu.value,
            self.photoionization_block_references,
        )
        stimulated_recombination_rate_coeff = pd.DataFrame(
            stimulated_recombination_rate_coeff,
            index=self.photoionization_index,
        )
        return stimulated_recombination_rate_coeff

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
            self.nu.value,
            self.photoionization_block_references,
        )
        photoionization_rate_coeff = pd.DataFrame(
            photoionization_rate_coeff,
            index=self.photoionization_index,
        )
        return photoionization_rate_coeff

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

        # Equation 15 Lucy 2003. Must be multiplied by factor Phi_ik from Eq 14
        stimulated_recombination_rate_coeff = (
            self.calculate_stimulated_recombination_rate_coeff(
                mean_intensity_photoionization_df,
                photoionization_boltzmann_factor,
            )
        )

        # Equation 16 Lucy 2003
        photoionization_rate_coeff = self.calculate_photoionization_rate_coeff(
            mean_intensity_photoionization_df,
        )

        return (
            photoionization_rate_coeff,
            stimulated_recombination_rate_coeff,
        )


class AnalyticCorrectedPhotoionizationCoeffSolver(
    SpontaneousRecombinationCoeffSolver
):
    def __init__(
        self,
        photoionization_cross_sections,
    ):
        super().__init__(photoionization_cross_sections)

    def calculate_mean_intensity_photoionization_df(
        self,
        radiation_field,
    ):
        """Calculates the mean intensity of the radiation field at each photoionization frequency.

        Parameters
        ----------
        radiation_field : RadiationField
            The radiation field.

        Returns
        -------
        pd.DataFrame
            DataFrame of mean intensities indexed by photoionization levels and
            columns of cells.
        """
        mean_intensity = radiation_field.calculate_mean_intensity(self.nu)
        return pd.DataFrame(
            mean_intensity,
            index=self.photoionization_cross_sections.index,
            columns=np.arange(len(radiation_field.temperature)),
        )

    def calculate_corrected_photoionization_rate_coeff(
        self,
        mean_intensity_photoionization_df,
        photoionization_boltzmann_factor,
        lte_level_population,
        level_population,
        lte_ion_population,
        ion_population,
    ):
        """
        Calculate the stimulated emission corrected photoionization rate coefficient.

        Parameters
        ----------
        mean_intensity_photoionization_df : pd.DataFrame
            A DataFrame of the mean intensity of the radiation field at each frequency

        Returns
        -------
        pd.DataFrame
            The calculated photoionization rate coefficient.

        Notes
        -----
        Equation 18 in Lucy 2003.
        """
        photoionization_rate_coeff = mean_intensity_photoionization_df.multiply(
            self.common_prefactor,
            axis=0,
        )

        # need to handle He and up. They have extra ionization states that
        # break the indexing.
        # Lucy 2003 Eq 18
        correction_factor = (
            1
            - (ion_population / lte_ion_population).values
            * (lte_level_population / level_population)
            * photoionization_boltzmann_factor
        )

        corrected_photoionization_rate_coeff = (
            photoionization_rate_coeff.multiply(correction_factor, axis=0)
        )

        corrected_photoionization_rate_coeff = integrate_array_by_blocks(
            corrected_photoionization_rate_coeff.values,
            self.nu.value,
            self.photoionization_block_references,
        )
        corrected_photoionization_rate_coeff = pd.DataFrame(
            corrected_photoionization_rate_coeff,
            index=self.photoionization_index,
        )
        return corrected_photoionization_rate_coeff

    def solve(
        self,
        dilute_blackbody_radiationfield_state,
        electron_temperature,
        lte_level_population,
        level_population,
        lte_ion_population,
        ion_population,
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
        photoionization_boltzmann_factor = pd.DataFrame(
            self.calculate_photoionization_boltzmann_factor(
                electron_temperature
            ),
            index=self.common_prefactor.index,
        )

        mean_intensity_photoionization_df = (
            self.calculate_mean_intensity_photoionization_df(
                dilute_blackbody_radiationfield_state
            )
        )
        # Equation 16 Lucy 2003
        corrected_photoionization_rate_coeff = (
            self.calculate_corrected_photoionization_rate_coeff(
                mean_intensity_photoionization_df,
                photoionization_boltzmann_factor,
                lte_level_population,
                level_population,
                lte_ion_population,
                ion_population,
            )
        )

        return corrected_photoionization_rate_coeff


class EstimatedPhotoionizationCoeffSolver:
    def __init__(
        self,
        level2continuum_edge_idx,
    ):
        self.level2continuum_edge_idx = level2continuum_edge_idx

    def solve(
        self,
        estimators_continuum,
        time_simulation,
        volume,
    ):
        """
        Solve for the continuum properties.

        Parameters
        ----------
        estimators_continuum : EstimatorsContinuum
            The Monte Carlo estimators for the continuum radiation field.
        time_simulation : float
            The simulation time.
        volume : float
            The volume of the cells.

        Returns
        -------
        ContinuumProperties
            The calculated continuum properties.

        Notes
        -----
        Lucy 2003 Eq 44, 45.
        """
        # TODO: the estimators are computed in the form epsilon_nu * distance * xsection / comoving_nu
        # with the stimulated recombination multiplied by a Boltzmann factor exp(-h * comoving_nu / k * electron_temp)
        # This is why this method does not match the one in AnalyticPhotoionizationCoeffSolver
        photoionization_normalization = (time_simulation * volume * H) ** -1

        photoionization_rate_coeff = bound_free_estimator_array2frame(
            estimators_continuum.photo_ion_estimator,
            self.level2continuum_edge_idx,
        )
        photoionization_rate_coeff *= photoionization_normalization

        stimulated_recombination_rate_coeff = bound_free_estimator_array2frame(
            estimators_continuum.stim_recomb_estimator,
            self.level2continuum_edge_idx,
        )
        stimulated_recombination_rate_coeff *= photoionization_normalization

        return photoionization_rate_coeff, stimulated_recombination_rate_coeff
