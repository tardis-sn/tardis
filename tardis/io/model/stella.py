import re
import dataclasses

import pandas as pd
from astropy import units as u


@dataclasses.dataclass
class STELLAModel:
    metadata: dict
    data: pd.DataFrame


HEADER_RE_STR = [
    ("\s+days post max Lbol\s+(\d+\.\d*)", "t_max"),
    ("\s+zones\s+(\d+)", "zones"),
    (
        "\s+inner boundary mass\s+(\d+\.\d+E[+-]\d+)\s+\d+\.\d+E[+-]\d+",
        "inner_boundary_mass",
    ),
    ("\s+total mass\s+(\d+\.\d+E[+-]\d+)\s+\d+\.\d+E[+-]\d+", "total_mass"),
]

DATA_START_ROW = 6

COLUMN_WITH_UNIT_RE = re.compile("(.+)\s+\((.+)\)")


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
            if line.strip().startswith("mass of cell"):
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
    data = pd.read_csv(
        fname,
        delim_whitespace=True,
        skiprows=DATA_START_ROW,
        header=None,
        index_col=0,
    )
    data.columns = column_names
    return STELLAModel(metadata, data)
