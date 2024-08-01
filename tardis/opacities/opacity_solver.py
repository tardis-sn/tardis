from tardis.opacities.tau_sobolev import calculate_sobolev_line_opacity
from tardis.opacities.opacity_state import OpacityState, opacity_state_initialize
import numpy as np

class OpacitySolver(object):

    def __init__(self, opacities_config):

        self.config = opacities_config
        self.line_interaction_type = opacities_config['line_interaction_type']
        self.disable_line_scattering = opacities_config['disable_line_scattering']

    def solve(self, legacy_plasma) -> OpacityState:

        tau_sobolev = calculate_sobolev_line_opacity(
            legacy_plasma.atomic_data.lines,
            legacy_plasma.level_number_density,
            legacy_plasma.time_explosion,
            legacy_plasma.stimulated_emission_factor,
        )

        if self.disable_line_scattering:
            tau_sobolev *= 0.0

        if self.line_interaction_type == "scatter":
            # to adhere to data types, we must have an array of minimum size 1
            array_size = 1
            transition_probabilities = np.zeros(
                (array_size, array_size), dtype=np.float64
            )  # to adhere to data types
            line2macro_level_upper = np.zeros(array_size, dtype=np.int64)
            macro_block_references = np.zeros(array_size, dtype=np.int64)
            transition_type = np.zeros(array_size, dtype=np.int64)
            destination_level_id = np.zeros(array_size, dtype=np.int64)
            transition_line_id = np.zeros(array_size, dtype=np.int64)

