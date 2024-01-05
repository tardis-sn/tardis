from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np
import pandas as pd
from astropy import units as u

import tardis.constants as const
from tardis.io.atom_data import AtomData
from tardis.model.radiation_field_state import (
    DiluteBlackBodyRadiationFieldState,
)
from tardis.montecarlo.estimators.util import (
    bound_free_estimator_array2frame,
    integrate_array_by_blocks,
)
from tardis.montecarlo.estimators.radfield_mc_estimators import (
    RadiationFieldMCEstimators,
)
from tardis.plasma.properties.continuum_processes import PhotoIonBoltzmannFactor

H = const.h.cgs.value


class EstimatedContinuumPropertiesABC(ABC):
    @abstractmethod
    def calculate_photo_ionization_rate_coefficient(self):
        pass

    @abstractmethod
    def calculate_stimulated_recomb_rate_coefficient(self):
        pass


class MCContinuumPropertiesSolver(EstimatedContinuumPropertiesABC):
    def __init__(
        self,
        atom_data: AtomData,
    ):
        self.atom_data = atom_data

    def solve(
        self,
        radfield_mc_estimators: RadiationFieldMCEstimators,
        time_simulation,
        volume,
    ):
        photo_ion_norm_factor = (time_simulation * volume * H) ** -1

        photo_ionization_rate_coefficient = bound_free_estimator_array2frame(
            radfield_mc_estimators.photo_ion_estimator,
            self.atom_data.level2continuum_edge_idx,
        )
        photo_ionization_rate_coefficient *= photo_ion_norm_factor

        stimulated_recomb_rate_factor = bound_free_estimator_array2frame(
            radfield_mc_estimators.stim_recomb_estimator,,
            self.atom_data.level2continuum_edge_idx,
        )
        stimulated_recomb_rate_factor *= photo_ion_norm_factor

        return ContinuumProperties(stimulated_recomb_rate_factor, photo_ionization_rate_coefficient)


class DiluteBlackBodyContinuumPropertiesSolver(EstimatedContinuumPropertiesABC):
    def __init__(self, atom_data: AtomData) -> None:
        self.atom_data = atom_data

    def solve(
        self,
        dilute_blackbody_radiationfield_state: DiluteBlackBodyRadiationFieldState,
        t_electrons: u.Quantity,
    ):
        photo_ion_boltzmann_factor = PhotoIonBoltzmannFactor.calculate(
            self.atom_data.photoionization_data, t_electrons
        )
        mean_intensity_photo_ion_df = (
            self.calculate_mean_intensity_photo_ion_table(
                dilute_blackbody_radiationfield_state
            )
        )

        photo_ion_rate_coeff = self.calculate_photo_ionization_rate_coefficient(
            mean_intensity_photo_ion_df
        )
        stimulated_recomb_rate_coeff = (
            self.calculate_stimulated_recomb_rate_coefficient(
                mean_intensity_photo_ion_df,
                photo_ion_boltzmann_factor,
            )
        )

        return ContinuumProperties(
            stimulated_recomb_rate_coeff, photo_ion_rate_coeff
        )

    def calculate_photo_ionization_rate_coefficient(
        self,
        mean_intensity_photo_ion_df: pd.DataFrame,
    ):
        gamma = mean_intensity_photo_ion_df.multiply(
            4.0
            * np.pi
            * self.atom_data.photoionization_data.x_sect
            / (self.atom_data.photoionization_data.nu * H),
            axis=0,
        )
        gamma = integrate_array_by_blocks(
            gamma.values,
            self.atom_data.photoionization_data.nu.values,
            self.atom_data.photo_ion_block_references,
        )
        gamma = pd.DataFrame(gamma, index=self.atom_data.photo_ion_unique_index)
        return gamma

    def calculate_stimulated_recomb_rate_coefficient(
        self,
        mean_intensity_photo_ion_df: pd.DataFrame,
        photo_ion_boltzmann_factor: np.ndarray,
    ):
        alpha_stim_factor = (
            mean_intensity_photo_ion_df * photo_ion_boltzmann_factor
        )
        alpha_stim_factor = alpha_stim_factor.multiply(
            4.0
            * np.pi
            * self.atom_data.photoionization_data.x_sect
            / (self.atom_data.photoionization_data.nu * H),
            axis=0,
        )
        alpha_stim_factor = integrate_array_by_blocks(
            alpha_stim_factor.values,
            self.atom_data.photoionization_data.nu.values,
            self.atom_data.photo_ion_block_references,
        )

        alpha_stim_factor = pd.DataFrame(
            alpha_stim_factor, index=self.atom_data.photo_ion_unique_index
        )

        return alpha_stim_factor

    def calculate_mean_intensity_photo_ion_table(
        self,
        dilute_blackbody_radiationfield_state: DiluteBlackBodyRadiationFieldState,
    ):
        mean_intensity = (
            dilute_blackbody_radiationfield_state.calculate_mean_intensity(
                self.atom_data.photoionization_data.nu.values
            )
        )
        mean_intensity_df = pd.DataFrame(
            mean_intensity,
            index=self.atom_data.photoionization_data.index,
            columns=np.arange(
                len(dilute_blackbody_radiationfield_state.t_radiative)
            ),
        )
        return mean_intensity_df


"""
class PhotoIonBoltzmannFactor(ProcessingPlasmaProperty):

    Attributes
    ----------
    boltzmann_factor_photo_ion : pandas.DataFrame, dtype float


    outputs = ("boltzmann_factor_photo_ion",)

    def calculate(self, photo_ion_cross_sections, t_electrons):
        nu = photo_ion_cross_sections["nu"].values

        boltzmann_factor = np.exp(-nu[np.newaxis].T / t_electrons * (H / K_B))
        return boltzmann_factor
"""


@dataclass
class ContinuumProperties:
    stimulated_recomb_rate_factor: pd.DataFrame
    photo_ionization_rate_coefficient: pd.DataFrame
