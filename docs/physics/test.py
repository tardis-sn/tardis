# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %% [markdown]
# ## Gamma-ray energy deposition

# %%
import numpy as np
import astropy.units as u

from tardis.energy_input.base import GXPacket, main_gamma_ray_loop
from tardis.energy_input.calculate_opacity import (
    compton_opacity_calculation,
    photoabsorption_opacity_calculation,
    pair_creation_opacity_calculation,
    kappa_calculation,
)
from tardis.energy_input.gamma_ray_grid import (
    distance_trace,
    move_gamma_ray,
    mass_distribution,
)
from tardis.energy_input.gamma_ray_interactions import scatter_type

# %% [markdown]
# Euler-Rodrigues rotation matrix for Compton Scattering rotation
# %% [markdown]
# ### Main loop
#
# Generates a simple 1D ejecta and a list of gamma-ray objects.
#
# Runs packets of gamma-rays through the ejecta. Handles interactions by calling the appropriate function.
#
# Appends deposited energy and output energy to 2 different lists. Currently no binning.

# %%
import matplotlib.pyplot as plt
from tardis.model import Radial1DModel
from tardis.io.config_reader import Configuration

num_packets = 500
num_packets = int(num_packets)

np.random.seed(1)

config = Configuration.from_yaml(
    "../../tardis/io/tests/data/tardis_configv1_density_exponential_nebular.yml"
)
config.model.structure.velocity.start = 1 * u.km / u.s
# config.model.structure.velocity.stop = 2e4 * u.km / u.s
# config.supernova.time_explosion = 240 * u.day
config.model.structure.density.rho_0 = 5e2 * u.g / (u.cm ** 3)
# config.model.structure.density.v_0 = 5e3 * u.km / u.s
# config.model.structure.density.time_0 = 240 * u.day

model = Radial1DModel.from_config(config)

(ejecta_energy_df, ejecta_plot_energy_df) = main_gamma_ray_loop(
    num_packets,
    model,
)

ejecta_energy = ejecta_energy_df["energy_input"]
ejecta_energy_r = ejecta_energy_df["energy_input_r"]
energy_input_time = ejecta_energy_df["energy_input_time"]
energy_input_type = ejecta_energy_df["energy_input_type"]


# %%
fig = plt.figure(dpi=150, facecolor="w")
ax = fig.add_subplot(111)

velocity = (
    np.array(ejecta_energy_r)
    / (config.supernova.time_explosion * u.day.to(u.s))
) * u.cm.to(u.km)
# ax.vlines(radii/np.max(model.r_outer.value), 0, 3500, linestyle="--", color="gray", alpha=0.2)
scatter = ax.scatter(
    np.array(ejecta_energy_r) / np.max(model.v_outer.value),
    np.array(ejecta_energy),
    c=energy_input_type,
    s=1,
    alpha=0.5,
)

ax.set_xlabel("r/R")
ax.set_ylabel("E (keV)")
ax.semilogy()


# %%
fig, ax = plt.subplots(
    subplot_kw={"projection": "polar"}, dpi=150, facecolor="w"
)
plot = ax.scatter(
    np.array(ejecta_energy_theta)[np.array(energy_input_type) == 1],
    np.array(ejecta_energy_r)[np.array(energy_input_type) == 1]
    / np.max(model.v_outer[:].value),
    c=np.array(ejecta_energy)[np.array(energy_input_type) == 1],
    s=1,
    alpha=0.5,
    cmap="plasma",
)
cbar = plt.colorbar(plot, ax=ax)


# %%
fig = plt.figure(dpi=150, facecolor="w")
ax = fig.add_subplot(111)
ax.hist(np.array(ejecta_energy_r) / np.max(model.v_outer[:].value), bins=200)
# ax.semilogy()
ax.set_xlim(0, 1)

# %% [markdown]
# # Total Mass

# %%
np.sum((model.density * model.volume)).to(u.M_sun)

# %% [markdown]
# # Density Profile

# %%
fig = plt.figure(dpi=150, facecolor="w")
plt.semilogy(model.r_middle / np.max(model.r_outer), model.density)
plt.plot(0, 0)


# %%
sum(escape_energy) / (sum(escape_energy) + sum(ejecta_energy))


# %%


# %%
