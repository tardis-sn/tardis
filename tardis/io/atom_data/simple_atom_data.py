from dataclasses import dataclass

import pandas as pd


@dataclass
class SimpleAtomData:
    data: pd.DataFrame
