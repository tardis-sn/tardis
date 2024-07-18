import re
import pandas as pd
from astropy import units as u
from pathlib import Path
import dataclasses


@dataclasses.dataclass
class CMFGENModel:
    metadata: dict
    data: pd.DataFrame


HEADER_RE_STR = [
    (r"t0:\s+(\d+\.\d+)+\s+day", "t0"),
]

COLUMN_ROW = 1
UNIT_ROW = 2
DATA_START_ROW = 3


def read_cmfgen_model(fname):
    """
    Read in a CMFGEN model file and return the data and model

    Parameters
    ----------

    fname : str

    Returns
    -------
    model : CMFGENModel

    """
    header_re = [re.compile(re_str[0]) for re_str in HEADER_RE_STR]
    metadata = {}
    with open(fname) as fh:
        for i, line in enumerate(fh):
            if i < len(HEADER_RE_STR):
                header_re_match = header_re[i].match(line)
                metadata[HEADER_RE_STR[i][1]] = header_re_match.group(1)
            elif i == COLUMN_ROW:
                if "Index" in line:
                    column_names = re.split(r"\s", line.strip())
                    column_names = [
                        col.lower().replace(" ", "_") for col in column_names
                    ]
                    column_names = column_names[
                        1:
                    ]  # Remove Index from column names
                else:
                    raise ValueError(
                        '"Index" is required in the Cmfgen input file to infer columns'
                    )
            elif i == UNIT_ROW:
                units = re.split(r"\s", line.strip())
                units = units[1:]  # Remove index column
                for col, unit in zip(column_names, units):
                    if u.Unit(unit) == "":  # dimensionless
                        continue
                    metadata[f"{col}_unit"] = u.Unit(unit)
                break

        metadata["t0"] = float(metadata["t0"]) * u.day
    data = pd.read_csv(
        fname,
        delim_whitespace=True,
        skiprows=DATA_START_ROW,
        header=None,
        index_col=0,
    )
    data.columns = column_names
    return CMFGENModel(metadata, data)
