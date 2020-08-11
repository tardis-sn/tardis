import os

from matplotlib import colors

import tardis.util.base
from tardis.io.atom_data import base
from matplotlib import pyplot as plt

import numpy as np
import pandas as pd

# Making 2 Figures for ionization balance and level populations

plt.figure(1).clf()
ax1 = plt.figure(1).add_subplot(111)

plt.figure(2).clf()
ax2 = plt.figure(2).add_subplot(111)

# expanding the tilde to the users directory
atom_fname = os.path.expanduser("~/.tardis/si_kurucz.h5")

# reading in the HDF5 File
atom_data = base.AtomData.from_hdf(atom_fname)

# The atom_data needs to be prepared to create indices. The Class needs to know which atomic numbers are needed for the
# calculation and what line interaction is needed (for "downbranch" and "macroatom" the code creates special tables)
atom_data.prepare_atom_data([14], "scatter")

# Initializing the NebularPlasma class using the from_abundance class method.
# This classmethod is normally only needed to test individual plasma classes
# Usually the plasma class just gets the number densities from the model class
nebular_plasma = plasma.NebularPlasma.from_abundance(
    10000, 0.5, {"Si": 1}, 1e-13, atom_data, 10.0
)


# Initializing a dataframe to store the ion populations  and level populations for the different temperatures
ion_number_densities = pd.DataFrame(index=nebular_plasma.ion_populations.index)
level_populations = pd.DataFrame(
    index=nebular_plasma.level_populations.loc[14, 1].index
)
t_rads = np.linspace(2000, 20000, 100)

# Calculating the different ion populations and level populuatios for the given temperatures
for t_rad in t_rads:
    nebular_plasma.update_radiationfield(t_rad, w=1.0)
    # getting total si number density
    si_number_density = nebular_plasma.number_density.get_value(14)
    # Normalizing the ion populations
    ion_density = nebular_plasma.ion_populations / si_number_density
    ion_number_densities[t_rad] = ion_density

    # normalizing the level_populations for Si II
    current_level_population = (
        nebular_plasma.level_populations.loc[14, 1]
        / nebular_plasma.ion_populations.loc[14, 1]
    )
    # normalizing with statistical weight
    current_level_population /= atom_data.levels.loc[14, 1].g

    level_populations[t_rad] = current_level_population

ion_colors = ["b", "g", "r", "k"]

for ion_number in [0, 1, 2, 3]:
    current_ion_density = ion_number_densities.loc[14, ion_number]
    ax1.plot(
        current_ion_density.index,
        current_ion_density.values,
        "%s-" % ion_colors[ion_number],
        label="Si %s W=1.0"
        % tardis.util.base.int_to_roman(ion_number + 1).upper(),
    )


# only plotting every 5th radiation temperature
t_rad_normalizer = colors.Normalize(vmin=2000, vmax=20000)
t_rad_color_map = plt.cm.ScalarMappable(norm=t_rad_normalizer, cmap=plt.cm.jet)

for t_rad in t_rads[::5]:
    ax2.plot(
        level_populations[t_rad].index,
        level_populations[t_rad].values,
        color=t_rad_color_map.to_rgba(t_rad),
    )
    ax2.semilogy()

# Calculating the different ion populations for the given temperatures with W=0.5
ion_number_densities = pd.DataFrame(index=nebular_plasma.ion_populations.index)
for t_rad in t_rads:
    nebular_plasma.update_radiationfield(t_rad, w=0.5)
    # getting total si number density
    si_number_density = nebular_plasma.number_density.get_value(14)
    # Normalizing the ion populations
    ion_density = nebular_plasma.ion_populations / si_number_density
    ion_number_densities[t_rad] = ion_density

    # normalizing the level_populations for Si II
    current_level_population = (
        nebular_plasma.level_populations.loc[14, 1]
        / nebular_plasma.ion_populations.loc[14, 1]
    )
    # normalizing with statistical weight
    current_level_population /= atom_data.levels.loc[14, 1].g

    level_populations[t_rad] = current_level_population

# Plotting the ion fractions

for ion_number in [0, 1, 2, 3]:
    print("w=0.5")
    current_ion_density = ion_number_densities.loc[14, ion_number]
    ax1.plot(
        current_ion_density.index,
        current_ion_density.values,
        "%s--" % ion_colors[ion_number],
        label="Si %s W=0.5"
        % tardis.util.base.int_to_roman(ion_number + 1).upper(),
    )

for t_rad in t_rads[::5]:
    ax2.plot(
        level_populations[t_rad].index,
        level_populations[t_rad].values,
        color=t_rad_color_map.to_rgba(t_rad),
        linestyle="--",
    )
    ax2.semilogy()

t_rad_color_map.set_array(t_rads)
cb = plt.figure(2).colorbar(t_rad_color_map)

ax1.set_xlabel("T [K]")
ax1.set_ylabel("Number Density Fraction")
ax1.legend()

ax2.set_xlabel("Level Number for Si II")
ax2.set_ylabel("Number Density Fraction")
cb.set_label("T [K]")

plt.show()
