from dataclasses import dataclass

import pandas as pd


@dataclass
class MacroAtomData:
    transition_probability_data: pd.DataFrame
    block_reference_data: pd.DataFrame
