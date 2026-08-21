from typing import Self

import numpy as np
import numpy.typing as npt
import pandas as pd

from tardis import constants as const
from tardis.opacities.continuum.continuum_state_numba import (
    ContinuumOpacityStateNumba,
)
from tardis.plasma.array_util import cumulative_integrate_array_by_blocks
from tardis.plasma.equilibrium.rates.heating_cooling_rates import (
    FreeFreeThermalRates,
)


class ContinuumOpacityState:
    """Store continuum data required for opacity computation."""

    def __init__(
        self,
        nu_i: pd.Series,
        level2continuum_idx: pd.Series,
        p_fb_deactivation: pd.DataFrame,
        photo_ion_cross_sections: pd.DataFrame,
        chi_bf: pd.DataFrame,
        ff_cooling_factor: npt.NDArray[np.float64],
        fb_emission_cdf: pd.DataFrame,
        photo_ion_idx: pd.DataFrame,
    ) -> None:
        """
        Store the continuum state required for opacity computation.

        Parameters
        ----------
        nu_i : pd.DataFrame
            frequencies for the bound-free thresholds
        level2continuum_idx : pd.DataFrame
            mapping from levels to the continuum
        p_fb_deactivation : pd.DataFrame
            probabilities of free-bound deactivation channels
        photo_ion_cross_sections : pd.DataFrame
            Photoionization cross sections
        chi_bf : pd.DataFrame
            Bound-free opacities
        ff_cooling_factor : np.ndarray
            free-free cooling factor
        fb_emission_cdf : pd.DataFrame
            free-bound emission cumulative distribution function
        photo_ion_idx : pd.DataFrame
            photoionization indices
        """
        self.nu_i = nu_i
        self.level2continuum_idx = level2continuum_idx
        self.p_fb_deactivation = p_fb_deactivation
        self.photo_ion_cross_sections = photo_ion_cross_sections
        self._chi_bf = chi_bf
        self.ff_cooling_factor = ff_cooling_factor
        self.fb_emission_cdf = fb_emission_cdf
        self.photo_ion_idx = photo_ion_idx

    @classmethod
    def from_equilibrium(
        cls,
        nu_i: pd.Series,
        level2continuum_idx: pd.Series,
        photo_ion_cross_sections: pd.DataFrame,
        photo_ion_idx: pd.DataFrame,
        ion_number_density: pd.DataFrame,
        level_number_density: pd.DataFrame,
        electron_densities: pd.Series,
        t_electrons: npt.ArrayLike,
        level_to_continuum_saha_factor: pd.DataFrame,
        radiative_recombination_rate: pd.DataFrame,
    ) -> Self:
        """Build transport continuum data from an accepted equilibrium state.

        Parameters
        ----------
        nu_i : pandas.Series
            Bound-free threshold frequencies.
        level2continuum_idx : pandas.Series
            Mapping from bound levels to continuum indices.
        photo_ion_cross_sections : pandas.DataFrame
            Photoionization frequencies and cross sections.
        photo_ion_idx : pandas.DataFrame
            Photoionization macro-atom transition indices.
        ion_number_density : pandas.DataFrame
            Accepted ion populations.
        level_number_density : pandas.DataFrame
            Accepted level populations.
        electron_densities : pandas.Series
            Accepted electron densities.
        t_electrons : numpy.typing.ArrayLike
            Accepted electron temperatures.
        level_to_continuum_saha_factor : pandas.DataFrame
            LTE bound-level population per next-ion and electron population.
        radiative_recombination_rate : pandas.DataFrame
            Accepted radiative recombination coefficients.

        Returns
        -------
        ContinuumOpacityState
            Continuum data formatted for transport.
        """
        photo_ion_index = radiative_recombination_rate.index
        upper_ion_index = pd.MultiIndex.from_arrays(
            [
                photo_ion_index.get_level_values("atomic_number"),
                photo_ion_index.get_level_values("ion_number") + 1,
            ],
            names=["atomic_number", "ion_number"],
        )
        recombination_cooling = (
            radiative_recombination_rate
            * electron_densities.to_numpy()
            * ion_number_density.loc[upper_ion_index].to_numpy()
        )
        cooling_totals = recombination_cooling.sum(axis=0).replace(0.0, 1.0)
        p_fb_deactivation = recombination_cooling.divide(cooling_totals, axis=1)
        p_fb_deactivation.index = level2continuum_idx.loc[
            p_fb_deactivation.index
        ].to_numpy()
        p_fb_deactivation.index.name = "continuum_idx"
        p_fb_deactivation.sort_index(inplace=True)

        frequencies = photo_ion_cross_sections.nu.to_numpy()
        boltzmann_factor = np.exp(
            -frequencies[:, None]
            / np.asarray(t_electrons)[None, :]
            * (const.h.cgs.value / const.k_B.cgs.value)
        )
        frequency_upper_ion_index = pd.MultiIndex.from_arrays(
            [
                photo_ion_cross_sections.index.get_level_values(
                    "atomic_number"
                ),
                photo_ion_cross_sections.index.get_level_values("ion_number")
                + 1,
            ],
            names=["atomic_number", "ion_number"],
        )
        stimulated_recombination_population = (
            level_to_continuum_saha_factor.loc[
                photo_ion_cross_sections.index
            ].to_numpy()
            * ion_number_density.loc[frequency_upper_ion_index].to_numpy()
            * electron_densities.to_numpy()
        )
        chi_bf = (
            level_number_density.loc[photo_ion_cross_sections.index]
            - stimulated_recombination_population * boltzmann_factor
        ).multiply(photo_ion_cross_sections.x_sect.to_numpy(), axis=0)

        ff_cooling_factor = (
            FreeFreeThermalRates()
            .heating_factor(
                ion_number_density,
                electron_densities.to_numpy(),
            )
            .to_numpy()
        )
        emission_integrand = (
            frequencies[:, None] ** 3
            * photo_ion_cross_sections.x_sect.to_numpy()[:, None]
            * boltzmann_factor
        )
        block_references = np.pad(
            photo_ion_cross_sections.nu.groupby(level=[0, 1, 2], sort=False)
            .count()
            .to_numpy()
            .cumsum(),
            [1, 0],
        ).astype(np.int64)
        emission_cdf = pd.DataFrame(
            cumulative_integrate_array_by_blocks(
                emission_integrand, frequencies, block_references
            ),
            index=photo_ion_cross_sections.index,
        ).fillna(0.0)

        return cls(
            nu_i,
            level2continuum_idx,
            p_fb_deactivation,
            photo_ion_cross_sections,
            chi_bf,
            ff_cooling_factor,
            emission_cdf,
            photo_ion_idx,
        )

    @classmethod
    def from_legacy_plasma(cls, plasma: object) -> Self:
        """
        Generate a ContinuumOpacityState object from a tardis BasePlasma.

        Parameters
        ----------
        plasma : tardis.plasma.BasePlasma
            legacy base plasma

        Returns
        -------
        ContinuumOpacityState
        """
        nu_i = plasma.nu_i
        level2continuum_idx = plasma.level2continuum_idx
        p_fb_deactivation = plasma.p_fb_deactivation
        photo_ion_cross_sections = plasma.photo_ion_cross_sections
        chi_bf = plasma.chi_bf
        ff_cooling_factor = plasma.ff_cooling_factor
        fb_emission_cdf = plasma.fb_emission_cdf
        photo_ion_idx = plasma.photo_ion_idx
        return cls(
            nu_i,
            level2continuum_idx,
            p_fb_deactivation,
            photo_ion_cross_sections,
            chi_bf,
            ff_cooling_factor,
            fb_emission_cdf,
            photo_ion_idx,
        )

    def to_numba(
        self,
        t_electrons: npt.NDArray[np.float64],
        macro_atom_state: object,
    ) -> ContinuumOpacityStateNumba:
        """Convert continuum data to the arrays used by transport.

        Parameters
        ----------
        t_electrons : numpy.ndarray
            Electron temperatures in each shell [K].

        Returns
        -------
        ContinuumOpacityStateNumba
            Array-backed continuum data used by IIP transport.
        """
        continuum_arrays = (
            self.bf_threshold_list_nu.values,
            np.ascontiguousarray(
                self.p_fb_deactivation.values.copy(), dtype=np.float64
            ),
            self.photo_ion_nu_threshold_mins.values,
            self.photo_ion_nu_threshold_maxs.values,
            self.photo_ion_block_references,
            self.chi_bf.values,
            self.x_sect.values,
            self.phot_nus.values,
            (self.ff_cooling_factor / np.sqrt(t_electrons)).astype(
                np.float64
            ),
            self.emissivities.values,
        )
        if macro_atom_state.photo_ion_block_idx < 0:
            photo_ion_activation_idx = self.photo_ion_activation_idx.to_numpy(
                dtype=np.int64
            )
        else:
            photo_ion_activation_idx = np.full(
                len(self.level2continuum_idx),
                macro_atom_state.photo_ion_block_idx,
                dtype=np.int64,
            )
        return ContinuumOpacityStateNumba(
            *continuum_arrays,
            photo_ion_activation_idx,
            np.int64(macro_atom_state.k_packet_idx),
            macro_atom_state.absorbing_probability_matrix,
        )

    @property
    def bf_threshold_list_nu(self) -> pd.Series:
        """
        List of Bound-Free Threshold Frequencies

        Returns
        -------
        pd.DataFrame
        """
        return self.nu_i.loc[self.level2continuum_idx.index]

    @property
    def phot_nus(self) -> pd.Series:
        """
        Frequencies corresponding to Photoionization Cross Sections

        Returns
        -------
        pd.DataFrame
        """
        return self.photo_ion_cross_sections.nu.loc[
            self.level2continuum_idx.index
        ]

    @property
    def photo_ion_block_references(self) -> npt.NDArray[np.int64]:
        """Photoionization Block References

        Returns
        -------
        np.ndarray
        """
        return np.pad(
            self.phot_nus.groupby(level=[0, 1, 2], sort=False)
            .count()
            .values.cumsum(),
            [1, 0],
        )

    @property
    def photo_ion_nu_threshold_mins(self) -> pd.Series:
        """
        Minimum Edges of the photoionization threshold frequencies

        Returns
        -------
        pd.DataFrame
        """
        return self.phot_nus.groupby(level=[0, 1, 2], sort=False).first()

    @property
    def photo_ion_nu_threshold_maxs(self) -> pd.Series:
        """
        Maximum Edges of the photoionization threshold frequencies

        Returns
        -------
        pd.DataFrame
        """
        return self.phot_nus.groupby(level=[0, 1, 2], sort=False).last()

    @property
    def x_sect(self) -> pd.Series:
        """
        Photoionization Cross Sections mapped to the continuum indices

        Returns
        -------
        pd.DataFrame
        """
        return self.photo_ion_cross_sections.x_sect.loc[
            self.level2continuum_idx.index
        ]

    @property
    def chi_bf(self) -> pd.DataFrame:
        """
        Bound-Free Opacities indices corresponding to the continuum levels

        Returns
        -------
        pd.DataFrame
        """
        return self._chi_bf.loc[self.level2continuum_idx.index]

    @property
    def emissivities(self) -> pd.DataFrame:
        """
        Free-bound Emissivities corresponding to the continuum levels

        Returns
        -------
        pd.DataFrame
        """
        return self.fb_emission_cdf.loc[self.level2continuum_idx.index]

    @property
    def photo_ion_activation_idx(self) -> pd.Series:
        """
        Index corresponding to photoionization activation

        Returns
        -------
        pd.DataFrame
        """
        return self.photo_ion_idx.loc[
            self.level2continuum_idx.index, "destination_level_idx"
        ]
