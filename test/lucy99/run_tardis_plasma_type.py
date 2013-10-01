#testing simple packet propagation with tardis and the normal kurucz dataset


import os
import pylab
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from tardis import atomic, model, montecarlo_multizone, synspec, config_reader
import numpy as np




def plot_tardis():
    radial1d_mdl = model.Radial1DModel(config_object, atom_data)
    radial1d_mdl.create_packets()
    out_nu, out_energy, jsestimator, nubarestimator = montecarlo_multizone.montecarlo_radial1d(radial1d_mdl)
    spectrum = synspec.get_lambda_spec(out_nu, out_energy, 500*1e-8, 20000*1e-8, samples=1000)
    ax.plot(np.log10(spectrum.wave/1e4), np.log10(spectrum.flux), label=config_object.plasma_type)
    ax.set_xlim(-1,0)
    ax.set_xlabel(r'$\log10{\lambda} [\log10{\mu m}]$')
    ax.set_ylabel(r'$\log10{Flux}$')
    return out_nu, out_energy, jsestimator, nubarestimator

#checking if kurucz_atom.h5 exists
kurucz_h5 = os.path.expanduser('~/.tardis/kurucz_atom.h5')
if not os.path.exists(kurucz_h5):
    raise IOError('kurucz_atom.h5 not found in %s' % kurucz_h5)


fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

config_object = config_reader.read_config('tardis_config_scatter.ini')

config_object.plasma_type='lte'
config_object.single_run_packets = int(1e5)

atom_data = atomic.AtomData.from_hdf5(kurucz_h5, use_macro_atom=True, use_zeta_data=True)

plot_tardis()
config_object.plasma_type='nebular'
out_nu, out_energy, jsestimator, nubarestimator = plot_tardis()

ax.legend()
fig.savefig('tardis_scatter_plasma_types.pdf')
fig.savefig('tardis_scatter_plasma_types.png')









