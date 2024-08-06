from dataclasses import dataclass

import pandas as pd


@dataclass
class CollisionData:
    data: pd.DataFrame
    temperatures: pd.DataFrame
    yg: pd.DataFrame
