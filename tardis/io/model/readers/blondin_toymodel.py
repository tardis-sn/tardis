import re

import numpy as np
import pandas as pd
import yaml
from astropy import units as u

from tardis.util.base import parse_quantity

PATTERN_REMOVE_BRACKET = re.compile(r"\[.+\]")
T0_PATTERN = re.compile("tend = (.+)\n")


def read_blondin_toymodel(fname):
    """
    Reading the Blondin toy-model format and returns a dictionary and a
    dataframe

    Parameters
    ----------
    fname : str
        path or filename to blondin toymodel

    Returns
    -------
    blondin_dict : dict
        dictionary containing most of the meta data of the model
    blondin_csv : pandas.DataFrame
        DataFrame containing the csv part of the toymodel
    """
    with open(fname, "r") as fh:
        for line in fh:
            if line.startswith("#idx"):
                break
        else:
            raise ValueError(
                "File {0} does not conform to Toy Model format as it does "
                "not contain #idx"
            )
    columns = [
        PATTERN_REMOVE_BRACKET.sub("", item) for item in line[1:].split()
    ]

    raw_blondin_csv = pd.read_csv(
        fname,
        # The argument `delim_whitespace` was changed to `sep`
        #   because the first one is deprecated since version 2.2.0.
        #   The regular expression means: the separation is one or
        #   more spaces together (simple space, tabs, new lines).
        sep=r"\s+",
        comment="#",
        header=None,
        names=columns,
    )
    raw_blondin_csv.set_index("idx", inplace=True)

    blondin_csv = raw_blondin_csv.loc[
        :,
        [
            "vel",
            "dens",
            "temp",
            "X_56Ni0",
            "X_Ti",
            "X_Ca",
            "X_S",
            "X_Si",
            "X_O",
            "X_C",
        ],
    ]
    rename_col_dict = {
        "vel": "velocity",
        "dens": "density",
        "temp": "t_electron",
    }
    rename_col_dict.update({item: item[2:] for item in blondin_csv.columns[3:]})
    rename_col_dict["X_56Ni0"] = "Ni56"
    blondin_csv.rename(columns=rename_col_dict, inplace=True)
    blondin_csv.iloc[:, 3:] = blondin_csv.iloc[:, 3:].divide(
        blondin_csv.iloc[:, 3:].sum(axis=1), axis=0
    )

    # changing velocities to outer boundary
    new_velocities = 0.5 * (
        blondin_csv.velocity.iloc[:-1].values
        + blondin_csv.velocity.iloc[1:].values
    )
    new_velocities = np.hstack(
        (new_velocities, [2 * new_velocities[-1] - new_velocities[-2]])
    )
    blondin_csv["velocity"] = new_velocities

    with open(fname, "r") as fh:
        t0_string = T0_PATTERN.findall(fh.read())[0]

    t0 = parse_quantity(t0_string.replace("DAYS", "day"))
    blondin_dict = {}
    blondin_dict["model_density_time_0"] = str(t0)
    blondin_dict["description"] = f"Converted {fname} to csvy format"
    blondin_dict["tardis_model_config_version"] = "v1.0"
    blondin_dict_fields = [
        dict(
            name="velocity",
            unit="km/s",
            desc="velocities of shell outer bounderies.",
        )
    ]
    blondin_dict_fields.append(
        dict(name="density", unit="g/cm^3", desc="mean density of shell.")
    )
    blondin_dict_fields.append(
        dict(name="t_electron", unit="K", desc="electron temperature.")
    )

    for abund in blondin_csv.columns[3:]:
        blondin_dict_fields.append(
            dict(name=abund, desc=f"Fraction {abund} abundance")
        )
    blondin_dict["datatype"] = {"fields": blondin_dict_fields}

    return blondin_dict, blondin_csv


def convert_blondin_toymodel(
    in_fname, out_fname, v_inner, v_outer, conversion_t_electron_rad=None
):
    """
    Parameters
    ----------
    in_fname : str
        input toymodel file
    out_fname : str
        output csvy file
    conversion_t_electron_rad : float or None
        multiplicative conversion factor from t_electron to t_rad.
        if `None` t_rad is not calculated
    v_inner : float or astropy.unit.Quantity
        inner boundary velocity. If float will be interpreted as km/s
    v_outer : float or astropy.unit.Quantity
        outer boundary velocity. If float will be interpreted as km/s
    """
    blondin_dict, blondin_csv = read_blondin_toymodel(in_fname)
    blondin_dict["v_inner_boundary"] = str(u.Quantity(v_inner, u.km / u.s))
    blondin_dict["v_outer_boundary"] = str(u.Quantity(v_outer, u.km / u.s))

    if conversion_t_electron_rad is not None:
        blondin_dict["datatype"]["fields"].append(
            {
                "desc": "converted radiation temperature "
                f"using multiplicative factor={conversion_t_electron_rad}",
                "name": "t_rad",
                "unit": "K",
            }
        )

        blondin_csv["t_rad"] = (
            conversion_t_electron_rad * blondin_csv.t_electron
        )

    csvy_file = f"---\n{yaml.dump(blondin_dict, default_flow_style=False)}\n---\n{blondin_csv.to_csv(index=False)}"

    with open(out_fname, "w") as fh:
        fh.write(csvy_file)
