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
LYMAN_CONTINUUM_INDEX = (1, 0, 0)


def _suppress_lyman_continuum(
    rate_coefficients: pd.DataFrame,
) -> pd.DataFrame:
    """Apply the configured H I ground-edge suppression policy."""
    if LYMAN_CONTINUUM_INDEX in rate_coefficients.index:
        rate_coefficients.loc[LYMAN_CONTINUUM_INDEX] = 0.0
    return rate_coefficients


def apply_lte_level_to_ion_factor(
    raw_rate_coefficients: pd.DataFrame,
    lte_level_to_ion_factor: pd.DataFrame,
) -> pd.DataFrame:
    """Apply the Lucy detailed-balance factor once to a raw coefficient.

    Parameters
    ----------
    raw_rate_coefficients : pandas.DataFrame
        Density-independent recombination coefficient indexed by photoionizing
        levels, with columns corresponding to shells.
    lte_level_to_ion_factor : pandas.DataFrame
        Lucy phi_ik detailed-balance factor with the same level index and shell
        columns.

    Returns
    -------
    pandas.DataFrame
        Recombination coefficient after applying ``phi_ik`` exactly once.
    """
    aligned_lte_level_to_ion_factor = lte_level_to_ion_factor.reindex(
        raw_rate_coefficients.index
    )
    if aligned_lte_level_to_ion_factor.isna().any().any():
        raise ValueError(
            "lte_level_to_ion_factor does not cover all recombination levels"
        )
    return raw_rate_coefficients.multiply(
        aligned_lte_level_to_ion_factor, axis=0
    )


def calculate_corrected_photoionization_rate_coeff(
    photoionization_rate_coefficients: pd.DataFrame,
    stimulated_recombination_rate_coefficients: pd.DataFrame,
    lte_level_population: pd.DataFrame,
    level_population: pd.DataFrame,
    lte_ion_population: pd.DataFrame,
    ion_population: pd.DataFrame,
) -> pd.DataFrame:
    """Apply the legacy stimulated-emission correction to raw coefficients.

    Parameters
    ----------
    photoionization_rate_coefficients : pandas.DataFrame
        Raw density-independent photoionization coefficients.
    stimulated_recombination_rate_coefficients : pandas.DataFrame
        Raw density-independent stimulated-recombination coefficients.
    lte_level_population, level_population : pandas.DataFrame
        LTE and current level populations.
    lte_ion_population, ion_population : pandas.DataFrame
        LTE and current ion populations.

    Returns
    -------
    pandas.DataFrame
        Population-corrected photoionization coefficient.

    This is an explicit rate assembly operation. The two input coefficients
    remain independent, density-free quantities and are not modified.
    """
    photoionization_index = photoionization_rate_coefficients.index
    lte_level_population = lte_level_population.reindex(photoionization_index)
    level_population = level_population.reindex(photoionization_index)
    lte_ion_population = _align_photoionization_population(
        lte_ion_population, photoionization_index
    )
    ion_population = _align_photoionization_population(
        ion_population, photoionization_index
    )
    population_ratio = (ion_population / lte_ion_population) * (
        lte_level_population / level_population
    )
    return photoionization_rate_coefficients - (
        population_ratio * stimulated_recombination_rate_coefficients
    )


def _align_photoionization_population(
    population: pd.DataFrame,
    photoionization_index: pd.MultiIndex,
) -> pd.DataFrame:
    """Align a level- or ion-indexed population with photoionizing levels."""
    if population.index.nlevels == photoionization_index.nlevels:
        return population.reindex(photoionization_index)

    destination_ions = pd.MultiIndex.from_arrays(
        [
            photoionization_index.get_level_values("atomic_number"),
            photoionization_index.get_level_values("ion_number") + 1,
        ],
        names=["atomic_number", "ion_number"],
    )
    aligned = population.reindex(destination_ions)
    aligned.index = photoionization_index
    return aligned


class SpontaneousRecombinationCoeffSolver:
    """Calculate raw spontaneous-recombination coefficients."""

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
        """Return the common Lucy coefficient prefactor.

        Used to multiply with both spontaneous recombination and
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

        return _suppress_lyman_continuum(
            spontaneous_recombination_rate_coeff_df
        )


class AnalyticPhotoionizationCoeffSolver(SpontaneousRecombinationCoeffSolver):
    """Calculate raw analytic photoionization coefficients."""

    def __init__(
        self,
        photoionization_cross_sections,
    ):
        super().__init__(photoionization_cross_sections)

    def calculate_mean_intensity_photoionization_df(
        self,
        dilute_blackbody_radiationfield_state,
    ):
        """Calculate the mean intensity at each photoionization frequency.

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
        """Prepare grouped ionization and recombination coefficients.

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
            _suppress_lyman_continuum(photoionization_rate_coeff),
            _suppress_lyman_continuum(stimulated_recombination_rate_coeff),
        )


class AnalyticCorrectedPhotoionizationCoeffSolver(
    SpontaneousRecombinationCoeffSolver
):
    """Assemble the legacy population-corrected photoionization coefficient."""

    def __init__(
        self,
        photoionization_cross_sections,
    ):
        super().__init__(photoionization_cross_sections)

    def calculate_mean_intensity_photoionization_df(
        self,
        radiation_field,
    ):
        """Calculate the mean intensity at each photoionization frequency.

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
        photoionization_index = self.photoionization_index
        lte_level_population = lte_level_population.reindex(
            photoionization_index
        )
        level_population = level_population.reindex(photoionization_index)
        lte_ion_population = _align_photoionization_population(
            lte_ion_population, photoionization_index
        )
        ion_population = _align_photoionization_population(
            ion_population, photoionization_index
        )

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
        """Prepare grouped ionization and recombination coefficients.

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

        return _suppress_lyman_continuum(corrected_photoionization_rate_coeff)


class EstimatedPhotoionizationCoeffSolver:
    """Expose raw Monte Carlo photoionization coefficients."""

    def __init__(
        self,
        level2continuum_edge_idx,
    ):
        self.level2continuum_edge_idx = level2continuum_edge_idx

    def solve(
        self,
        estimators_continuum: dict[str, object],
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Return raw Monte Carlo photoionization coefficients.

        Parameters
        ----------
        estimators_continuum : EstimatorsContinuum
            The Monte Carlo estimators for the continuum radiation field.

        Returns
        -------
        tuple[pandas.DataFrame, pandas.DataFrame]
            Raw photoionization and stimulated-recombination coefficients.

        Notes
        -----
        Lucy 2003 Eq 44, 45.
        """
        # TODO: the estimators are computed in the form epsilon_nu * distance * xsection / comoving_nu
        # with the stimulated recombination multiplied by a Boltzmann factor exp(-h * comoving_nu / k * electron_temp)
        # This is why this method does not match the one in AnalyticPhotoionizationCoeffSolver
        # photoionization_normalization = (time_simulation * volume * H) ** -1

        photoionization_rate_coeff = bound_free_estimator_array2frame(
            estimators_continuum["photoionization_rate_estimator"],
            self.level2continuum_edge_idx,
        )

        stimulated_recombination_rate_coeff = bound_free_estimator_array2frame(
            estimators_continuum["stimulated_recombination_rate_estimator"],
            self.level2continuum_edge_idx,
        )

        return (
            _suppress_lyman_continuum(photoionization_rate_coeff),
            _suppress_lyman_continuum(stimulated_recombination_rate_coeff),
        )

    def solve_corrected(
        self,
        estimators_continuum: dict[str, object],
        level_population: pd.DataFrame,
        lte_level_population: pd.DataFrame,
        ion_population: pd.DataFrame,
        lte_ion_population: pd.DataFrame,
    ) -> pd.DataFrame:
        """Assemble the legacy corrected photoionization coefficient.

        Parameters
        ----------
        estimators_continuum : dict[str, object]
            Monte Carlo bound-free estimators.
        level_population, lte_level_population : pandas.DataFrame
            Current and LTE level populations.
        ion_population, lte_ion_population : pandas.DataFrame
            Current and LTE ion populations.

        Returns
        -------
        pandas.DataFrame
            Population-corrected photoionization coefficient.
        """
        photoionization_rate_coeff, stimulated_recombination_rate_coeff = (
            self.solve(estimators_continuum)
        )
        return calculate_corrected_photoionization_rate_coeff(
            photoionization_rate_coeff,
            stimulated_recombination_rate_coeff,
            lte_level_population,
            level_population,
            lte_ion_population,
            ion_population,
        )
