import re
from dataclasses import dataclass

import pandas as pd
from astropy import units as u


@dataclass
class STELLAModel:
    metadata: dict
    data: pd.DataFrame


HEADER_RE_STR = [
    (r"\s+days post max Lbol\s+(\d+\.\d*)", "t_max"),
    (r"\s+zones\s+(\d+)", "zones"),
    (
        r"\s+inner boundary mass\s+(\d+\.\d+E[+-]\d+)\s+\d+\.\d+E[+-]\d+",
        r"inner_boundary_mass",
    ),
    (r"\s+total mass\s+(\d+\.\d+E[+-]\d+)\s+\d+\.\d+E[+-]\d+", "total_mass"),
]

DATA_START_ROW = 5
COLUMN_WITH_UNIT_RE = re.compile(r"(.+)\s+\((.+)\)")


def read_stella_model(fname):
    """
    Read in a STELLA model file and return the data and model

    Parameters
    ----------

    fname : str

    Returns
    -------
    model : STELLAModel

    """
    header_re = [re.compile(re_str[0]) for re_str in HEADER_RE_STR]
    metadata = {}
    with open(fname) as fh:
        for i, line in enumerate(fh):
            if i < len(HEADER_RE_STR):
                header_re_match = header_re[i].match(line)
                metadata[HEADER_RE_STR[i][1]] = header_re_match.group(1)
            elif i == DATA_START_ROW:
                if "mass of cell" in line:
                    column_names_raw = re.split(r"\s{3,}", line.strip())
                    break
                else:
                    raise ValueError(
                        '"mass of cell" is required in the Stella input file to infer columns'
                    )

        metadata["t_max"] = float(metadata["t_max"]) * u.day
        metadata["zones"] = int(metadata["zones"])
        metadata["inner_boundary_mass"] = (
            float(metadata["inner_boundary_mass"]) * u.g
        )
        metadata["total_mass"] = float(metadata["total_mass"]) * u.g
    column_names = []
    for column_name in column_names_raw:
        if (column_match := COLUMN_WITH_UNIT_RE.match(column_name)) is not None:
            column_name, column_unit = column_match.groups()
            column_name = column_name.lower().replace(" ", "_")
            metadata[f"{column_name}_unit"] = u.Unit(column_unit)
        else:
            column_name = column_name.lower().replace(" ", "_")
        column_names.append(column_name)
    # +1 because there is a missing line between columns
    # and the actual data
    data = pd.read_csv(
        fname,
        # The argument `delim_whitespace` was changed to `sep`
        #   because the first one is deprecated since version 2.2.0.
        #   The regular expression means: the separation is one or
        #   more spaces together (simple space, tabs, new lines).
        sep=r"\s+",
        skiprows=DATA_START_ROW + 1,
        header=None,
        index_col=0,
    )
    data.columns = column_names
    return STELLAModel(metadata, data)
