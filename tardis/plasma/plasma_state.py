from tardis.io.util import HDFWriterMixin
from tardis.transport.montecarlo.configuration import montecarlo_globals
import numpy as np

# NOTE: we are missing the plasma atomic_data properties.  Are any of these ever updated?
# NOTE: Should we have individual states for different interactions?  Should plasma_state be a union or list of states?
# TODO: MacroAtomData should set macro_block_references even when not using continuum for consistency
# NOTE: Currently all of the properties returned are as they are in the plasma, so dataframes, arrays, units, etc is ambiguous
# NOTE: continuum ff_opacity factor should be computed by the opacity state


def get_default(cls, attr, default=None):

    return getattr(cls, attr) if hasattr(cls, attr) else default


class PlasmaState(HDFWriterMixin):
    """Current State of the Plasma Required for Opacity Computation"""

    hdf_name = "plasma_state"

    hdf_properties = [
        "electron_densities",
        "electron_temperature",
        "level_number_density",
        "stimulated_emission_factor",
        "transition_probabilities",
        "continuum_state",
        "macro_atom_state",
    ]

    def __init__(
        self,
        electron_densities,
        electron_temperature,
        level_number_density,
        stimulated_emission_factor,
        macro_atom_state=None,
        continuum_state=None,
    ):

        self.electron_densities = electron_densities
        self.electron_temperature = electron_temperature

        # These required for calculating sobolev opacities
        self.level_number_density = level_number_density
        self.stimulated_emission_factor = stimulated_emission_factor

        self.macro_atom_state = macro_atom_state
        self.continuum_state = continuum_state

    @classmethod
    def from_plasma(cls, plasma, continuum_state=None, macro_atom_state=None):

        electron_densities = plasma.electron_densities
        electron_temperature = plasma.t_electrons

        # These required for calculating sobolev opacities
        level_number_density = plasma.level_number_density
        stimulated_emission_factor = plasma.stimulated_emission_factor

        # If line interaction type is not scatter

        return cls(
            electron_densities,
            electron_temperature,
            level_number_density,
            stimulated_emission_factor,
            continuum_state,
            macro_atom_state,
        )


class MacroAtomState(HDFWriterMixin):

    """Current State of the MacroAtom"""

    hdf_name = "macro_atom_state"

    hdf_properties = [
        "transition_probabilities",
        "transition_type",
        "destination_level_id",
        "transition_line_id",
        "macro_block_references",
    ]

    def __init__(
        self,
        transition_probabilities,
        transition_type,
        destination_level_id,
        transition_line_id,
        macro_block_references,
    ):

        self.transition_probabilities = transition_probabilities
        self.transition_type = transition_type
        self.destination_level_id = destination_level_id
        self.transition_line_id = transition_line_id
        self.macro_block_references = macro_block_references

    @classmethod
    def from_plasma(cls, plasma):

        transition_probabilities = plasma.macro_atom_data["transition_type"]
        transition_type = plasma.macro_atom_data["transition_type"]
        destination_level_id = plasma.macro_atom_data["destination_level_idx"]
        transition_line_id = plasma.macro_atom_data["lines_idx"]

        if (
            montecarlo_globals.CONTINUUM_PROCESSES_ENABLED
        ):  # TODO: Unify this in the plasma solver
            macro_block_references = plasma.macro_block_references
        else:
            macro_block_references = plasma.atomic_data.macro_atom_references[
                "block_references"
            ]

        return cls(
            transition_probabilities,
            transition_type,
            destination_level_id,
            transition_line_id,
            macro_block_references,
        )


class ContinuumState(HDFWriterMixin):
    """Current State of the Continuum Required for Opacity Computation"""

    hdf_name = "continuum_state"

    hdf_properties = [
        "nu_i",
        "level2continuum_idx",
        "p_fb_deactivation",
        "photo_ion_cross_sections",
        "chi_bf",
        "ff_cooling_factor",
        "fb_emission_cdf",
        "photo_ion_idx",
        "k_packet_idx",
    ]

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
    def from_plasma(cls, plasma):

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
