from abc import ABC, abstractmethod

import numpy as np
import pandas as pd

import tardis.constants as const
from tardis.model.radiation_field_state import (
    DiluteBlackBodyRadiationFieldState,
)
from tardis.montecarlo.estimators.util import bound_free_estimator_array2frame

H = const.h.cgs.value


class EstimatedContinuumPropertiesABC(ABC):
    @abstractmethod
    def calculate_photo_ionization_rate_coefficient(self):
        pass

    @abstractmethod
    def calculate_stimulated_recomb_rate_coefficient(self):
        pass


class EstimatedMCContinuumProperties(EstimatedContinuumPropertiesABC):
    def __init__(
        self, estimator_statistics, level2continuum_idx, time_simulation, volume
    ):
        self.estimator_statistics = estimator_statistics
        self.level2continuum_idx = level2continuum_idx

        # initializing the lazy properties
        self._gamma_estimator_df = None
        self._alpha_stimulated_estimator_df = None
        self.photo_ion_norm_factor = (time_simulation * volume * H) ** -1

    def calculate_photo_ionization_rate_coefficient(self):
        return self.gamma_estimator_df * self.photo_ion_norm_factor.value

    def calculate_stimulated_recomb_rate_coefficient(self):
        alpha_stim = (
            self.alpha_stimulated_estimator_df * self.photo_ion_norm_factor
        )

        # does that need the new Phis or the old ones?
        alpha_stim *= phi_ik.loc[alpha_stim.index]
        return alpha_stim

    @property
    def alpha_stimulated_estimator_df(self):
        if self._alpha_stimulated_estimator_df is None:
            self._alpha_stimulated_estimator_df = (
                bound_free_estimator_array2frame(
                    self.estimator_statistics.stim_recomb_estimator,
                    self.level2continuum_idx,
                )
            )
        return self._alpha_stimulated_estimator_df

    @property
    def gamma_estimator_df(self):
        if self._gamma_estimator_df is None:
            self._gamma_estimator_df = bound_free_estimator_array2frame(
                self.estimator_statistics.photo_ion_estimator,
                self.level2continuum_idx,
            )
        return self._gamma_estimator_df


class EstimatedDiluteBlackBodyContinuumProperties(
    EstimatedContinuumPropertiesABC
):
    def __init__(
        self,
        dilute_blackbody_radiationfield_state: DiluteBlackBodyRadiationFieldState,
    ) -> None:
        self.dilute_blackbody_radiationfield_state = (
            dilute_blackbody_radiationfield_state
        )

    def calculate_photo_ionization_rate_coefficient(
        self,
        photo_ion_cross_sections,
        photo_ion_block_references,
        photo_ion_index,
    ):
        nu = photo_ion_cross_sections["nu"]
        cross_section = photo_ion_cross_sections["x_sect"]
        mean_intensity = (
            self.dilute_blackbody_radiationfield_state.calculate_mean_intensity(
                photo_ion_cross_sections.nu
            )
        )
        mean_intensity_df = pd.DataFrame(
            mean_intensity,
            index=photo_ion_cross_sections.index,
            columns=np.arange(len(self.dilute_blackbody_radiationfield_state)),
        )

        gamma = j_nus.multiply(4.0 * np.pi * cross_section / nu / H, axis=0)
        gamma = integrate_array_by_blocks(
            gamma.values, nu.values, photo_ion_block_references
        )
        gamma = pd.DataFrame(gamma, index=photo_ion_index)
        return gamma
