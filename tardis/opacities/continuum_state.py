import numpy as np


class ContinuumState:
    """Current State of the Continuum Required for Opacity Computation"""

    def __init__(
        self,
        nu_i,
        level2continuum_idx,
        p_fb_deactivation,
        photo_ion_cross_sections,
        chi_bf,
        ff_cooling_factor,
        fb_emission_cdf,
        photo_ion_idx,
        k_packet_idx,
    ):

        self.nu_i = nu_i
        self.level2continuum_idx = level2continuum_idx
        self.p_fb_deactivation = p_fb_deactivation
        self.photo_ion_cross_sections = photo_ion_cross_sections
        self._chi_bf = chi_bf
        self.ff_cooling_factor = ff_cooling_factor
        self.fb_emission_cdf = fb_emission_cdf
        self.photo_ion_idx = photo_ion_idx
        self.k_packet_idx = k_packet_idx

    @classmethod
    def from_legacy_plasma(cls, plasma):

        nu_i = plasma.nu_i
        level2continuum_idx = plasma.level2continuum_idx
        p_fb_deactivation = plasma.p_fb_deactivation
        photo_ion_cross_sections = plasma.photo_ion_cross_sections
        chi_bf = plasma.chi_bf
        ff_cooling_factor = plasma.ff_cooling_factor
        fb_emission_cdf = plasma.fb_emission_cdf
        photo_ion_idx = plasma.photo_ion_idx
        k_packet_idx = plasma.k_packet_idx

        return cls(
            nu_i,
            level2continuum_idx,
            p_fb_deactivation,
            photo_ion_cross_sections,
            chi_bf,
            ff_cooling_factor,
            fb_emission_cdf,
            photo_ion_idx,
            k_packet_idx,
        )

    @property
    def bf_threshold_list_nu(self):
        return self.nu_i.loc[self.level2continuum_idx.index]


    @property
    def phot_nus(self):
        return self.photo_ion_cross_sections.nu.loc[
            self.level2continuum_idx.index
        ]

    @property
    def photo_ion_block_references(self):

        return np.pad(
            self.phot_nus.groupby(level=[0, 1, 2], sort=False)
            .count()
            .values.cumsum(),
            [1, 0],
        )

    @property
    def photo_ion_nu_threshold_mins(self):

        return self.phot_nus.groupby(level=[0, 1, 2], sort=False).first()

    @property
    def photo_ion_nu_threshold_maxs(self):

        return self.phot_nus.groupby(level=[0, 1, 2], sort=False).last()

    @property
    def x_sect(self):
        return self.photo_ion_cross_sections.x_sect.loc[
            self.level2continuum_idx.index
        ]

    @property
    def chi_bf(self):
        return self._chi_bf.loc[self.level2continuum_idx.index]
 
    @property
    def emissivities(self):
        return self.fb_emission_cdf.loc[self.level2continuum_idx.index]

    @property
    def photo_ion_activation_idx(self):
        return self.photo_ion_idx.loc[
            self.level2continuum_idx.index, "destination_level_idx"
        ]