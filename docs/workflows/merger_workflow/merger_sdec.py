import os
import pdb
import shutil
from io import StringIO
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import constants as const
from astropy import units as u
from plot_manga import plot_manga_profiles

from tardis.io.configuration.config_reader import Configuration
from tardis.visualization import SDECPlotter
from tardis.workflows.nonhomologous_tardis_workflow import (
    NonhomologousTARDISWorkflow,
)
from tardis.workflows.util import get_tau_integ


def load_csvy_dataframe(path):
    path = Path(path)
    lines = path.read_text().splitlines()
    if not lines or lines[0].strip() != "---":
        return pd.read_csv(path)

    csv_start = next(
        idx + 1
        for idx, line in enumerate(lines[1:], start=1)
        if line.strip() == "---"
    )
    return pd.read_csv(StringIO("\n".join(lines[csv_start:])))



# Load in the MANGA data
time_explosion = 100 # days
fname = f"bins_{time_explosion}days.dat"
nangles = 19
nradii = 599
angle_bin = 9

index = pd.MultiIndex.from_product([range(nangles), range(nradii)], names=['angle bin', 'radius bin'])
data = pd.DataFrame(np.loadtxt(fname, comments=['#', '\n']), columns=["r", "dens", "temp", "v_r"], index=index)

# Change selected_bin to highlight a different line of sight
fig = plot_manga_profiles(data, nangles, selected_bin=angle_bin, day=time_explosion)
axes = fig.axes

# Set label for the MANGA data
axes[5].lines[0].set_label('MANGA 1D data for highlighted bin')

# Load the processed 1D profile from the csvy
merger_df = load_csvy_dataframe("manga_merger.csvy")

merger_df["velocity_cms"] = (
    merger_df["velocity"].to_numpy() * u.km / u.s
).to_value(u.cm / u.s)
merger_df["radius_cm"] = (
    merger_df["velocity"].to_numpy() * u.km / u.s * time_explosion * u.day
).to_value(u.cm)

# Overplot processed csvy data
axes[3].plot(merger_df["radius_cm"], merger_df["velocity_cms"], c="C1", lw=5, alpha=0.3)
axes[4].plot(merger_df["radius_cm"], merger_df["density"], c="C1", lw=5, alpha=0.3)
axes[5].plot(merger_df["radius_cm"], merger_df["t_rad"], c="C1", lw=5, alpha=0.3, label='CSVY data used for TARDIS sim')

# Plot fiducial tardis example profiles for comparison
config = Configuration.from_yaml('../../tardis_example.yml')
workflow = NonhomologousTARDISWorkflow(config)
v_fiducial = np.linspace(config.model.structure.velocity.start, config.model.structure.velocity.stop, config.model.structure.velocity.num)
dens_fiducial = workflow.simulation_state.composition.density
temp_fiducial = workflow.simulation_state.t_radiative
r_fiducial = (v_fiducial*time_explosion*u.day).to(u.cm)
axes[3].plot(r_fiducial, v_fiducial.to(u.cm/u.s), c='C3')
axes[4].plot(r_fiducial, dens_fiducial, c='C3')
axes[5].plot(r_fiducial, temp_fiducial, c='C3', label='tardis_example.yml (for reference)')
axes[3].set_xscale('log')
axes[4].set_xscale('log')
axes[5].set_xscale('log')
axes[5].legend(frameon=False)



# Create workflow from config
config_yaml_path = Path("manga_merger_config.yml")
config = Configuration.from_yaml(str(config_yaml_path))

workflow = NonhomologousTARDISWorkflow(config, csvy=True)


# Calculate optical depth by integrating in from the surface
opacity_states = workflow.solve_opacity()
tau_integ = get_tau_integ(
    workflow.plasma_solver,
    opacity_states['opacity_state'],
    workflow.simulation_state,
)
tau_integ['rosseland'].to(1)

# Electron scattering only: kappa_es = n_e * sigma_T
n_e = workflow.plasma_solver.electron_densities.values * u.cm**-3
kappa_es = (n_e * const.sigma_T).to(u.cm**-1)

dr = (workflow.simulation_state.geometry.r_outer - workflow.simulation_state.geometry.r_inner).to(u.cm)
dtau_es = (kappa_es * dr).to(u.dimensionless_unscaled)

# Integrate from outside in (cumsum reversed)
tau_es_integrated = np.cumsum(dtau_es[::-1])[::-1]

photosphere_radius_idx = np.argwhere(tau_es_integrated <= 2./3.)[0]
try:
    r_phot = merger_df['radius_cm'][1:][photosphere_radius_idx].values[0]
    tau_inner = tau_es_integrated[0].value[0]
except:
    r_phot = merger_df['radius_cm'][1]
    tau_inner = tau_es_integrated[0].value


print(f"photopshere at {r_phot:.2e}")
print(f"tau at photosphere is {tau_inner:.3f}")

lum_logLsun = np.log10(config.supernova.luminosity_requested.to(u.Lsun).value)
writepath = Path(f"./100day_nonhomologous_{lum_logLsun:.2f}-logLsun_{tau_inner:.3f}-tau-inner/")
os.makedirs(writepath.absolute(), exist_ok=True)

plt.savefig(writepath/"manga_profile.png")
plt.close()

shutil.copy("manga_merger.csvy", writepath)
shutil.copy(config_yaml_path, writepath)

fig = plt.figure()
plt.plot(merger_df['radius_cm'][1:], tau_integ['rosseland'].to(1), label='Rosseland Mean')
plt.plot(merger_df['radius_cm'][1:], tau_es_integrated, label='Electron Scattering')
plt.axvline(r_phot, c='k', ls=':')

plt.ylabel(r"$\tau$")
plt.xlabel("r (cm)")
plt.axhline(2./3, c='r', ls='--', label=r'$\tau = 2/3$')
plt.loglog()
plt.legend(frameon=False)
plt.savefig(writepath/"tau_photosphere.png")
plt.close()


# Now run the workflow and create an SDEC plot
workflow.converged = True
workflow.run()

fig = plt.figure()
plotter = SDECPlotter.from_workflow(workflow)
plotter.generate_plot_mpl(packets_mode="real", nelements=30, packet_wvl_range=[6400, 6700]*u.AA, show_modeled_spectrum=False)
plt.savefig(writepath/"sdec_Halpha.png")
plt.close()

fig = plt.figure()
plotter.generate_plot_mpl(packets_mode="real", nelements=30, packet_wvl_range=[2000, 20000]*u.AA, show_modeled_spectrum=False)
plt.savefig(writepath/"sdec_full.png")
plt.close()

fig = plt.figure()
plotter.generate_plot_mpl(packets_mode="real", nelements=30, packet_wvl_range=[2000, 20000]*u.AA, show_modeled_spectrum=False)

import pickle
with open(writepath/"sdec.p", "wb") as f:
    pickle.dump(fig, f)

pdb.set_trace()