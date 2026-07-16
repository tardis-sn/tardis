from dataclasses import dataclass

import numpy as np
import pandas as pd

from tardis import constants as const
from tardis.plasma.equilibrium.continuum_state import (
    EquilibriumContinuumState,
    _cumulative_integrate_by_blocks,
)


@dataclass
class EquilibriumContinuumOpacityState:
    """Continuum opacity data derived from equilibrium plasma quantities."""

    nu_i: pd.Series
    level2continuum_idx: pd.Series
    p_fb_deactivation: pd.DataFrame
    photo_ion_cross_sections: pd.DataFrame
    _chi_bf: pd.DataFrame
    ff_cooling_factor: np.ndarray
    fb_emission_cdf: pd.DataFrame
    photo_ion_idx: pd.DataFrame
    k_packet_idx: int

    @classmethod
    def from_plasma(
        cls, plasma, continuum_state: EquilibriumContinuumState
    ) -> EquilibriumContinuumOpacityState:
        """Build transport opacity data without legacy continuum properties."""
        photo_data = plasma.photo_ion_cross_sections
        photo_index = plasma.photo_ion_index
        level2continuum_idx = plasma.level2continuum_idx
        nu_i = plasma.nu_i

        upper_ion_index = pd.MultiIndex.from_arrays(
            [
                photo_index.get_level_values("atomic_number"),
                photo_index.get_level_values("ion_number") + 1,
            ],
            names=["atomic_number", "ion_number"],
        )
        recombination_cooling = (
            continuum_state.radiative_recombination_rate
            * plasma.electron_densities.values
            * plasma.ion_number_density.loc[upper_ion_index].values
        )
        cooling_totals = recombination_cooling.sum(axis=0).replace(0.0, 1.0)
        p_fb_deactivation = recombination_cooling.divide(cooling_totals, axis=1)
        continuum_indices = level2continuum_idx.loc[
            p_fb_deactivation.index
        ].values
        p_fb_deactivation = p_fb_deactivation.set_index(continuum_indices)
        p_fb_deactivation.index.name = "continuum_idx"
        p_fb_deactivation = p_fb_deactivation.sort_index()

        temperatures = plasma.t_electrons
        frequencies = photo_data.nu.to_numpy()
        boltzmann_factor = np.exp(
            -frequencies[:, np.newaxis]
            / temperatures[np.newaxis, :]
            * (const.h.cgs.value / const.k_B.cgs.value)
        )
        photo_levels = plasma.level_number_density.loc[photo_data.index]
        lte_photo_levels = plasma.lte_level_number_density.loc[photo_data.index]
        level_ratio = lte_photo_levels.divide(photo_levels)
        upper_ion_index_by_frequency = pd.MultiIndex.from_arrays(
            [
                photo_data.index.get_level_values("atomic_number"),
                photo_data.index.get_level_values("ion_number") + 1,
            ],
            names=["atomic_number", "ion_number"],
        )
        ion_ratio = plasma.ion_number_density.loc[
            upper_ion_index_by_frequency
        ].divide(
            plasma.lte_ion_number_density.loc[upper_ion_index_by_frequency]
        )
        chi_bf = (
            photo_levels
            * (1.0 - level_ratio.multiply(ion_ratio.values) * boltzmann_factor)
        ).multiply(photo_data.x_sect.to_numpy(), axis=0)

        charges = plasma.ion_number_density.index.get_level_values(
            "ion_number"
        ).to_numpy()
        ff_cooling_factor = (
            plasma.electron_densities
            * plasma.ion_number_density.multiply(charges**2, axis=0).sum()
        ).to_numpy()

        alpha_sp_energy = (
            frequencies[:, np.newaxis] ** 3
            * photo_data.x_sect.to_numpy()[:, np.newaxis]
            * boltzmann_factor
        )
        block_references = np.pad(
            photo_data.nu.groupby(level=[0, 1, 2], sort=False)
            .count()
            .to_numpy()
            .cumsum(),
            [1, 0],
        )
        fb_emission_cdf = pd.DataFrame(
            _cumulative_integrate_by_blocks(
                alpha_sp_energy, frequencies, block_references
            ),
            index=photo_data.index,
        )

        return cls(
            nu_i,
            level2continuum_idx,
            p_fb_deactivation,
            photo_data,
            chi_bf,
            ff_cooling_factor,
            fb_emission_cdf,
            plasma.photo_ion_idx,
            int(plasma.atomic_data.macro_atom_references.references_idx.max()),
        )

    @property
    def bf_threshold_list_nu(self) -> pd.Series:
        """Bound-free threshold frequencies in transport continuum order."""
        return self.nu_i.loc[self.level2continuum_idx.index]

    @property
    def phot_nus(self) -> pd.Series:
        """Photoionization frequencies in transport continuum order."""
        return self.photo_ion_cross_sections.nu.loc[
            self.level2continuum_idx.index
        ]

    @property
    def photo_ion_block_references(self) -> np.ndarray:
        """Block boundaries for the frequency-resolved photoionization data."""
        return np.pad(
            self.phot_nus.groupby(level=[0, 1, 2], sort=False)
            .count()
            .to_numpy()
            .cumsum(),
            [1, 0],
        )

    @property
    def photo_ion_nu_threshold_mins(self) -> pd.Series:
        """Lower frequency bound of each active photoionization block."""
        return self.phot_nus.groupby(level=[0, 1, 2], sort=False).first()

    @property
    def photo_ion_nu_threshold_maxs(self) -> pd.Series:
        """Upper frequency bound of each active photoionization block."""
        return self.phot_nus.groupby(level=[0, 1, 2], sort=False).last()

    @property
    def chi_bf_for_transport(self) -> pd.DataFrame:
        """Bound-free opacity contributions in transport continuum order."""
        return self._chi_bf.loc[self.level2continuum_idx.index]

    @property
    def x_sect(self) -> pd.Series:
        """Photoionization cross sections in transport continuum order."""
        return self.photo_ion_cross_sections.x_sect.loc[
            self.level2continuum_idx.index
        ]

    @property
    def emissivities(self) -> pd.DataFrame:
        """Free-bound emission CDFs in transport continuum order."""
        return self.fb_emission_cdf.loc[self.level2continuum_idx.index]

    @property
    def photo_ion_activation_idx(self) -> pd.Series:
        """Macro-atom activation indices for active continuum edges."""
        return self.photo_ion_idx.loc[
            self.level2continuum_idx.index, "destination_level_idx"
        ]

    @property
    def chi_bf(self) -> pd.DataFrame:
        """Alias for the transport-ordered bound-free opacity table."""
        return self.chi_bf_for_transport
