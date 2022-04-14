'''This is a code to realize the first project of M2D visualization, running a simulation with logging enabled with only one packet.'''
from tardis import run_tardis
from tardis.io.config_reader import Configuration
from tardis.io.atom_data.util import download_atom_data

# Downloading atom data
download_atom_data('kurucz_cd23_chianti_H_He')

# Reading the Configuration stored in `tardis_config_packet_tracking.yml` into␣ ↪config
config = Configuration.from_yaml("tardis_example.yml")

# Enabled with only one packet and allowing to track
config["montecarlo"]["no_of_packets"]=config["montecarlo"]["last_no_of_packets"]=1
config["montecarlo"]["no_of_virtual_packets"]=0 
config["montecarlo"]["iterations"]=2 
config["montecarlo"]["tracking"]["track_rpacket"]=True

# Running a simulation
sim = run_tardis(config, show_convergence_plots=False, show_progress_bars=False)

# Plot mu vs. radius
from matplotlib import pyplot as plt

fig, ax = plt.subplots()
ax.set_xlabel('r [cm]')
ax.set_ylabel('$\\rm \\mu$')
ax.plot(sim.runner.rpacket_tracker[0].r, sim.runner.rpacket_tracker[0].mu) 
plt.show()


