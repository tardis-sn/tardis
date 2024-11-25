from dataclasses import dataclass

import pandas as pd


@dataclass
class MacroAtomData:
    """Wrapper for macro atom related dataframes

    transition_probability_data : (pandas.DataFrame, np.array)
        index : numerical index
        columns : atomic_number, ion_number, source_level_number, destination_level_number,
        transition_line_id, transition_type, transition_probability;

    block_reference_data : (pandas.DataFrame, np.array)
        index : numerical index
        columns : atomic_number, ion_number, source_level_number, count_down, count_up, count_total.
        Refer to the docs: http://tardis.readthedocs.io/en/latest/physics/plasma/macroatom.html
    """

    transition_probability_data: pd.DataFrame
    block_reference_data: pd.DataFrame
