


class MacroAtomState:

    def __init__(self,
                 transition_probabilities,
                 macro_block_references,
                 transition_type,
                 destination_level_id,
                 transition_line_id,
                 ):

        self.transition_probabilities = transition_probabilities
        self.macro_block_references = macro_block_references
        self.transition_type = transition_type
        self.destination_level_id = destination_level_id
        self.transition_line_id = transition_line_id


    @classmethod
    def from_macro_atom_state(cls, macro_atom_state):

        return cls(    
            macro_atom_state.transition_probabilities,        
            macro_atom_state.macro_block_references,
            macro_atom_state.transition_type,
            macro_atom_state.destination_level_id,
            macro_atom_state.transition_line_id,
            )