# %%
from tardis import run_tardis
from tardis.io.config_reader import Configuration
import numpy as np
from astropy import units as u
import pandas as pd

config = Configuration.from_yaml("docs/models/examples/tardis_example.yml")
del config["config_dirname"]
config.model.abundances["O"] = 1.0
config.model.abundances["Mg"] = 0.0
config.model.abundances["Si"] = 0.0
config.model.abundances["Ar"] = 0.0
config.model.abundances["S"] = 0.0
config.model.abundances["Ca"] = 0.0
config.model.abundances["Fe"] = 0.0

config.montecarlo.iterations = 1

sim = run_tardis(config)

import tardis.energy_input.spencer_fano as sf
from collections import namedtuple


def read_colliondata(
    collionfilename="/home/afullard/spencer-fano/data/collion.txt",
):

    collionrow = namedtuple(
        "collionrow",
        [
            "atomic_number",
            "nelec",
            "n",
            "l",
            "ion_potential",
            "A",
            "B",
            "C",
            "D",
        ],
    )

    with open(collionfilename, "r") as collionfile:
        print(f"Collionfile: expecting {collionfile.readline().strip()} rows")
        dfcollion = pd.read_csv(
            collionfile,
            delim_whitespace=True,
            header=None,
            names=collionrow._fields,
        )
    dfcollion.eval("ion_number = atomic_number - nelec", inplace=True)

    return dfcollion


sim.plasma.ion_collision_data = read_colliondata()

# %%
(
    electron_spectrum,
    transitions_dict,
    energy_grid,
    ions,
    ion_populations,
    energy_deposition_density,
    ion_collision_data,
    electron_number_density,
    number_density,
) = sf.setup_solution(
    0.1,
    3000,
    5,
    temperature=10000,
    plasma=sim.plasma,
    energy_deposition_density=100,
)

# TODO get atomic data from TARDIS for transitions and levels


# %%
import matplotlib.pyplot as plt

plt.figure(dpi=150)
plt.semilogy(energy_grid.grid, electron_spectrum)
plt.ylim(0, 1e8)
plt.xlabel("Energy (eV)")
plt.ylabel("y(E)")


# %%
electron_spectrum


# %%
import tardis.energy_input.non_thermal_fractions as ntf

ntf.total_fractions(
    energy_grid,
    electron_spectrum,
    transitions_dict,
    ion_collision_data,
    ion_populations,
    energy_deposition_density,
    electron_number_density,
    number_density,
)


# %%
import matplotlib.pyplot as plt

level = x.index.get_level_values("level_number")
energy = x.energy.values

plt.plot(level, energy, ".")


# %%
sim.plasma.level_boltzmann_factor


# %%
sim.plasma.partition_function


# %%
sim.plasma.partition_function.loc[:, 0]


# %%
ion.name


# %%
atomic_transition_data = sim.plasma.atomic_data.macro_atom_data


# %%
atomic_number = 8
ion_number = 0

ion = atomic_transition_data.query(
    "atomic_number == @atomic_number and ion_number == @ion_number"
)


# %%
levels = sim.plasma.atomic_data.levels.query(
    "ion_number == 0" + " and atomic_number == 8"
)

ion = sim.plasma.atomic_data.macro_atom_data.query(
    "atomic_number == 8 and ion_number == 0"
)

print(len(ion))

src_energy = []
dest_energy = []

for level_id in ion["destination_level_number"]:
    dest_energy.append(levels.query("level_number == @level_id").energy.values)

for level_id in ion["source_level_number"]:
    src_energy.append(levels.query("level_number == @level_id").energy.values)


# %%
np.array(dest_energy) - np.array(src_energy) > 0


# %%
