from tardis.montecarlo.estimators.util import bound_free_estimator_array2frame
import tardis.constants as const

H = const.h.cgs.value


class ContinuumEstimators:
    def __init__(
        self, estimator_statistics, level2continuum_idx, time_simulation, volume
    ):
        self.estimator_statistics = estimator_statistics
        self.level2continuum_idx = level2continuum_idx

        # initializing the lazy properties
        self._gamma_estimator_df = None
        self._alpha_stimulated_estimator_df = None
        self.photo_ion_norm_factor = (time_simulation * volume * H) ** -1

    def calculate_gamma(self):
        return self.gamma_estimator_df * self.photo_ion_norm_factor.value

    def calculate_stimulated_recomb_rate_coeff(self):
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
