#testing simple packet propagation with tardis and the normal kurucz dataset


import os
import pylab
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from tardis import atomic, model_radial_oned, montecarlo_multizone, synspec
import numpy as np

kurucz_h5 = os.path.expanduser('~/.tardis/kurucz_atom.h5')
if not os.path.exists(kurucz_h5):
    raise IOError('kurucz_atom.h5 not found in %s' % kurucz_h5)

#checking if kurucz_atom.h5 exists

atom_data = atomic.AtomData.from_hdf5(kurucz_h5, use_macro_atom=True, use_zeta_data=True)
radial1d_mdl = model_radial_oned.Radial1DModel.from_config_file('tardis_config_scatter.ini', atom_data)
radial1d_mdl.create_packets()
out_nu, out_energy, jsestimator, nubarestimator = montecarlo_multizone.montecarlo_radial1d(radial1d_mdl)
spectrum = synspec.get_lambda_spec(out_nu, out_energy, 500*1e-8, 20000*1e-8, samples=1000)
fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.plot(np.log10(spectrum.wave/1e4), np.log10(spectrum.flux))
ax.set_xlim(-1,0)
ax.set_xlabel(r'$\log10{\lambda} [\log10{\mu m}]$')
ax.set_ylabel(r'$\log10{Flux}$')




