import re
import pandas as pd
from astropy import units as u
from pathlib import Path
import dataclasses


@dataclasses.dataclass
class STELLAModel:
    header: dict
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
    header = {}
    with open(fname) as fh:
        for i, line in enumerate(fh):
            header_re_match = header_re[i].match(line)
            header[HEADER_RE_STR[i][1]] = header_re_match.group(1)
            if i == len(header_re) - 1:
                break
        header["t_max"] = float(header["t_max"]) * u.day
        header["zones"] = int(header["zones"])
        header["inner_boundary_mass"] = (
            float(header["inner_boundary_mass"]) * u.g
        )
        header["total_mass"] = float(header["total_mass"]) * u.g

    data = pd.read_csv(fname, delim_whitespace=True, skiprows=7)
    return STELLAModel(header, data)
