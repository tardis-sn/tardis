from pathlib import Path  # noqa: I001

from tardis.io.atom_data import AtomData

import numpy as np
import pandas as pd
from astropy import units as u


def convert_asplund_to_mass_fraction(atom_data, abundances):

    # Only use species for which atom data actually exists
    species_requested = atom_data.levels.index.unique(level=0)

    # TEMPORARY: manually specify which species to include
    #include = [1, 2, 6, 7, 8]
    #species_requested = species_requested[species_requested.isin(include)]

    # Get atomic mass and abundance for each species
    masses = atom_data.atom_data.mass.loc[species_requested].values
    filtered_abundances = abundances.loc[abundances['atomic_number'].isin(species_requested)].reset_index(drop=True)

    # Convert to linear number density
    nX_over_nH = 10.0 ** (filtered_abundances.Value.values - 12.0)

    # Rescale to a total mass fraction instead of relative to hydrogen
    norm = (nX_over_nH * masses).sum()
    mass_fraction = nX_over_nH * masses / norm

    # Put data into a dataframe matching the format of the original abundance data
    mass_fraction_df = filtered_abundances.drop(columns="Value").copy()
    mass_fraction_df["MassFraction"] = mass_fraction

    return mass_fraction_df


def manga_to_csvy(data, abundances=None):
    v_km_s = (data["v_r"].values * u.cm / u.s).to(u.km / u.s).value# * 10 # TEMP: scaling factor to test larger velocities
    r_cm = data["r"].values
    rho = data["dens"].values
    t_rad = data["temp"].values

    # Generate entries in header for each elemental abundance
    abund_field_template = lambda x: f"    - name: {x}"

    if abundances is None:
        # Presume pure H
        abund_str = [abund_field_template("H")]
        abund_data = [1.0]
    else:
        abund_str = [
            abund_field_template(abundances.Element[i])
            for i in range(len(abundances))
        ]
        abund_data = [
            abundances.MassFraction[i] for i in range(len(abundances))
        ]

    csvy_yaml_header = f"""\
---
name: manga_merger_{time_explosion}day_bin{angle_bin}
model_density_time_0: {time_explosion} day
model_isotope_time_0: nan s
description: MANGA merger non-homologous model at {time_explosion} days, angle bin {angle_bin}
tardis_model_config_version: v1.0
datatype:
  fields:
    - name: velocity
      unit: km/s
      desc: velocities of shell boundaries (inner-to-outer)
    - name: radius
      unit: cm
      desc: actual radial positions of shell boundaries (non-homologous)
    - name: density
      unit: g/cm^3
      desc: shell density
    - name: t_rad
      unit: K
      desc: initial radiative temperature
{"\n".join(abund_str)}
v_inner_boundary: {v_km_s[0]:.6f} km/s
v_outer_boundary: {v_km_s[-1]:.6f} km/s
---
"""

    csvy_data = {
        "velocity": v_km_s,
        "radius": r_cm,
        "density": rho,
        "t_rad": t_rad,
    }
    csvy_data = pd.DataFrame(data=csvy_data)
    for species, frac in zip(abundances.Element, abundances.MassFraction):
        csvy_data[species] = frac

    csvy_content = csvy_yaml_header + "\n" + csvy_data.to_csv(index=False)
    csvy_path = Path("manga_merger.csvy")
    csvy_path.write_text(csvy_content)


# Load in the MANGA data
time_explosion = 100  # days
fname = f"bins_{time_explosion}days.dat"
nangles = 19
nradii = 599
angle_bin = 9

index = pd.MultiIndex.from_product(
    [range(nangles), range(nradii)], names=["angle bin", "radius bin"]
)
data = pd.DataFrame(
    np.loadtxt(fname, comments=["#", "\n"]),
    columns=["r", "dens", "temp", "v_r"],
    index=index,
)

data_1d = data.loc[angle_bin]
data_1d_cut = data_1d.loc[
    np.logical_and(data_1d["r"] > 1.57e14, data_1d["r"] < 3e14)
]
nshells = len(data_1d_cut) - 1

# Floor the temperature at 1e3 K
#data_1d_cut.loc[data_1d_cut["temp"] < 1e3, "temp"] = 1000.01
# TEMP: isothermal at 5000 K
data_1d_cut.loc[:, "temp"] = 5000.0

# TEMP: Create homologous data to compare with
data_homologous = data_1d_cut.copy(deep=True)
data_homologous['v_r'] = (data_homologous['r'] / (time_explosion * u.day.to(u.s)))
data_1d_cut['v_r'] = data_homologous['v_r']

# Load atom data and solar abundances
atom_data = AtomData.from_hdf("kurucz_cd23_cmfgen_H_He.h5")
asplund = pd.read_csv("asplund_2020_processed.csv")
abundances = convert_asplund_to_mass_fraction(atom_data, asplund)

manga_to_csvy(data_1d_cut, abundances)
