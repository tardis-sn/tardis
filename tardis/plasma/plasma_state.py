from tardis.io.util import HDFWriterMixin
from tardis.transport.montecarlo.configuration import montecarlo_globals
import numpy as np

# NOTE: we are missing the plasma atomic_data properties.  Are any of these ever updated?
# NOTE: Should we have individual states for different interactions?  Should plasma_state be a union or list of states?
# TODO: MacroAtomData should set macro_block_references even when not using continuum for consistency
# NOTE: Currently all of the properties returned are as they are in the plasma, so dataframes, arrays, units, etc is ambiguous
# NOTE: continuum ff_opacity factor should be computed by the opacity state


class PlasmaState(HDFWriterMixin):
    """Current State of the Plasma Required for Opacity Computation"""

    hdf_name = "plasma_state"

    hdf_properties = [
        "electron_densities",
        "electron_temperature",
        "level_number_density",
        "stimulated_emission_factor",
        "transition_probabilities",
    ]

    def __init__(
        self,
        electron_densities,
        electron_temperature,
        level_number_density,
        stimulated_emission_factor,
        macro_atom_state=None,
    ):

        self.electron_densities = electron_densities
        self.electron_temperature = electron_temperature

        # These required for calculating sobolev opacities
        self.level_number_density = level_number_density
        self.stimulated_emission_factor = stimulated_emission_factor

        self.macro_atom_state = macro_atom_state

    @classmethod
    def from_legacy_plasma(cls, plasma):

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
    def from_legacy_plasma(cls, plasma):

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
