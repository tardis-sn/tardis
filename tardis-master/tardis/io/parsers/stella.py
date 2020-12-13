import re
import pandas as pd
from astropy import units as u
import numpy as np


def read_stella_data(filename):
    with open(filename) as fh:
        col = fh.readlines()[5]
    col_names = re.split(r"\s{3,}", col.strip())
    col_names = [
        re.sub(r"\s\(.+\)", "", col_name).replace(" ", "_")
        for col_name in col_names
    ]
    data = pd.read_csv(
        filename, skiprows=7, delim_whitespace=True, names=col_names
    )

    # drop last row of data
    data = data.iloc[0:-1]
    return data
