import sys
import logging
import os

# Create a matplotlib backend that doesn't require a display
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from tardis.io.atom_data import download_atom_data
from tardis import run_tardis

logging.basicConfig(level=logging.INFO)

print("Downloading atom data...")
download_atom_data('kurucz_cd23_chianti_H_He_latest')

print("Running TARDIS simulation...")
# Run the simulation using the example config
# Make sure to give it the correct path to the yaml file
yaml_path = os.path.join("docs", "tardis_example.yml")
sim = run_tardis(yaml_path, 
                 virtual_packet_logging=True,
                 show_convergence_plots=False,
                 export_convergence_plots=False,
                 log_level="INFO")

print("Generating spectrum plot...")
# Save the spectrum plot
spectrum = sim.spectrum_solver.spectrum_real_packets
spectrum_virtual = sim.spectrum_solver.spectrum_virtual_packets
spectrum_integrated = sim.spectrum_solver.spectrum_integrated

plt.figure(figsize=(10, 6.5))
spectrum.plot(label="Normal packets")
spectrum_virtual.plot(label="Virtual packets")
spectrum_integrated.plot(label='Formal integral')
plt.xlim(500, 9000)
plt.title("TARDIS example model spectrum")
plt.xlabel(r"Wavelength [$\AA$]")
plt.ylabel(r"Luminosity density [erg/s/$\AA$]")
plt.legend()
plt.savefig("spectrum_output.png")

print("Simulation complete. Output saved to spectrum_output.png")
