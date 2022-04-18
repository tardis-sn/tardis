'''This is a code to realize the first project of M2D visualization, running a simulation with logging enabled with only one packet.'''
from tardis import run_tardis
from tardis.io.config_reader import Configuration
from tardis.io.atom_data.util import download_atom_data

# Reading the Configuration stored in `tardis_config_packet_tracking.yml` into␣ ↪config
config = Configuration.from_yaml("docs/tardis_example.yml")

# Enabled with only one packet and allowing to track
config["montecarlo"]["no_of_packets"]=config["montecarlo"]["last_no_of_packets"]=1
config["montecarlo"]["no_of_virtual_packets"]=0 
config["montecarlo"]["iterations"]=1 
config["montecarlo"]["tracking"]["track_rpacket"]=True

# Running a simulation
sim = run_tardis(config, show_convergence_plots=False, show_progress_bars=False)

# Plot mu vs. radius

import numpy as np
import matplotlib.pyplot as plt

# Using colormap to represent frequency
cm = plt.cm.get_cmap('gist_rainbow')
vmin = np.min(sim.runner.input_nu)*0.8
vmax = np.max(sim.runner.input_nu)*1.2

fig, ax = plt.subplots(figsize=[10,8])
ax.set_xlabel('r [cm]')
ax.set_ylabel('$\\rm \\mu$')
# Plotting trajectory
ax.plot(sim.runner.rpacket_tracker[0].r, sim.runner.rpacket_tracker[0].mu, c='gray') 
plt.scatter(sim.runner.rpacket_tracker[0].r, sim.runner.rpacket_tracker[0].mu, c=sim.runner.rpacket_tracker[0].nu, vmin=vmin, vmax=vmax, cmap=cm)
shell_length = np.size(sim.model.r_outer)
max_mu = np.max(sim.runner.rpacket_tracker[0].mu)
r_outer = sim.model.r_outer.to('cm').value
r_inner = sim.model.r_inner.to('cm').value
# Plotting shell boundaries
for i in range(shell_length):
    ax.axvline(r_outer[i], c='black', ls='--')
    ax.text(r_inner[i], max_mu, "%d"%(i+1))
plt.colorbar(label = '$\\nu [Hz]$')
plt.show()


