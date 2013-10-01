#testing simple packet propagation with tardis and the normal kurucz dataset


import os
import pylab
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from astropy import units

from tardis import atomic, model, montecarlo_multizone, synspec, config_reader
import numpy as np




def plot_tardis():
    radial1d_mdl = model.Radial1DModel(config_object, atom_data)
    radial1d_mdl.create_packets()
    out_nu, out_energy, jsestimator, nubarestimator = montecarlo_multizone.montecarlo_radial1d(radial1d_mdl)
    spectrum = synspec.get_lambda_spec(out_nu, out_energy, 500*1e-8, 20000*1e-8, samples=1000)
    ax.plot(np.log10(spectrum.wave/1e4), np.log10(spectrum.flux), label='emitted')
    ax.set_xlim(-1,0)
    ax.set_xlabel(r'$\log10{\lambda} [\log10{\mu m}]$')
    ax.set_ylabel(r'$\log10{Flux}$')
    return out_nu, out_energy, jsestimator, nubarestimator, radial1d_mdl

#checking if kurucz_atom.h5 exists
kurucz_h5 = os.path.expanduser('~/.tardis/kurucz_atom.h5')
if not os.path.exists(kurucz_h5):
    raise IOError('kurucz_atom.h5 not found in %s' % kurucz_h5)


fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

config_object = config_reader.read_config('tardis_config_scatter.ini')

config_object.plasma_type='nebular'
config_object.line_interaction_type='downbranch'
config_object.single_run_packets = int(1e5)

atom_data = atomic.AtomData.from_hdf5(kurucz_h5, use_macro_atom=True, use_zeta_data=True)

out_nu, out_energy, jsestimator, nubarestimator, radial1d_mdl = plot_tardis()

emission_filter = out_nu >= 0
input = synspec.get_lambda_spec(radial1d_mdl.packet_src.packet_nus, radial1d_mdl.packet_src.packet_energies, 500*1e-8, 20000*1e-8, samples=1000)

input_lambda = units.Unit('Hz').to('Angstrom', radial1d_mdl.packet_src.packet_nus, units.spectral())
output_lambda = units.Unit('Hz').to('Angstrom', out_nu, units.spectral())
fig.clf()
ax = fig.add_subplot(111)
ax.plot(input_lambda[emission_filter], output_lambda[emission_filter], 'b,')
ax.set_xlim(2000, 10000)
ax.set_ylim(2000, 10000)
ax.set_xlabel('Input Packet Wavelength')
ax.set_ylabel('Output Packet Wavelength')
ax.set_title('Input vs Output (Nebular Plasma & Downbranch)')

fig.savefig('tardis_nebular_downbranch_in_vs_out.pdf')
fig.savefig('tardis_nebular_downbranch_in_vs_out.png')









