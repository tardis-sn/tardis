import re
import pandas as pd
from astropy import units as u
import numpy as np
import yaml
from scipy import interpolate
import os


STELLA_CONFIG_TEMPLATE_FNAME = os.path.join(
    os.path.dirname(__file__), "tardis_stella_template.yml.j2"
)

STELLA_COL_MAPPER = {
    "mass of cell (g)": "mass_cell",
    "cell center m (g)": "mass_total",
    "cell center R (cm)": "r_center",
    "cell center v (cm/s)": "v_center",
    "avg density": "density",
    "radiation pressure": "pressure_rad",
    "avg temperature": "t_gas",
    "radiation temperature": "t_radiation",
    "avg opacity": "opacity",
    "outer edge m (g)": "mass_outer",
    "outer edge r (cm)": "r_outer",
}


def read_stella(fname):
    """
    Read a STELLA model into a pandas dataframe
    """
    stella_meta = {}
    with open(fname) as fh:
        for line in fh:
            if line.strip().startswith("days post max Lbol"):
                stella_meta["time_explosion"] = float(line[20:].strip()) * u.day
            if line.strip().startswith("mass"):
                break
    raw_columns = re.split("\s{3,}", line.strip())

    stella_model = pd.read_csv(
        fname, delim_whitespace=True, skiprows=6, names=raw_columns, index_col=0
    )
    stella_model.rename(columns=STELLA_COL_MAPPER, inplace=True)
    return stella_meta, stella_model


def convert_stella_to_csvy(fname, out_fname=None):
    """
    Convert STELLA to TARDIS file
    """
    stella_meta, stella_model = read_stella(fname)
    v_interpolator = interpolate.interp1d(
        stella_model["r_center"], stella_model["v_center"]
    )
    # velocity needed for future check against homology
    #velocity = v_interpolator(stella_model["r_outer"].values) * u.cm / u.s
    homologous_velocity = (
        stella_model.r_outer.values * u.cm / (stella_meta["time_explosion"])
    )
    csvy_meta = {
        "name": "STELLA transformation",
        "model_density_time_0": str(stella_meta["time_explosion"]),
        "model_isotope_time_0": str(stella_meta["time_explosion"]),
        "description": "parsed stella file",
        "tardis_model_config_version": "v1.0",
        "datatype": {"fields": []},
    }
    csvy_table = pd.DataFrame()
    csvy_table["density"] = stella_model.density
    csvy_meta["datatype"]["fields"].append(
        {
            "name": "density",
            "unit": "g/cm^3",
            "desc": "Average density from stella",
        }
    )

    csvy_table["velocity"] = homologous_velocity.to(u.km / u.s).value
    csvy_meta["datatype"]["fields"].append(
        {
            "name": "velocity",
            "unit": "km/s",
            "desc": "WARNING - THIS IS NOT THE STELLA VELOCITY but a homologous approximation given radius and t_explosion",
        }
    )

    abundances = stella_model.iloc[:, 12:-3]
    for isotope in abundances.columns:
        csvy_meta["datatype"]["fields"].append(
            {"name": isotope, "desc": "{0} mass fraction".format(isotope)}
        )

    csvy_table = csvy_table.join(abundances)
    if out_fname is None:
        out_fname = f"{fname}_stella2tardis.csvy"

    with open(out_fname, "w") as fh:
        fh.write('---\n')
        yaml.dump(csvy_meta, fh)
        fh.write('---\n')
        csvy_table.to_csv(fh, index=False)
