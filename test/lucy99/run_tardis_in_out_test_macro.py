#testing simple packet propagation with tardis and the normal kurucz dataset

######### IMPORTANT THIS SCRIPT WILL ONLY RUN WHEN WE ADD LINE_ID_IN AND LINE_ID_OUT TO MONTECARLO_MULTIZONE.pyx
####### this is a pure debug statement

import os
import pylab
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from astropy import units

from tardis import atomic, model_radial_oned, montecarlo_multizone, synspec, config_reader
import numpy as np




def plot_tardis():
    radial1d_mdl = model_radial_oned.Radial1DModel(config_object, atom_data)
    radial1d_mdl.create_packets()
    out_nu, out_energy, jsestimator, nubarestimator, line_id_in, line_id_out = montecarlo_multizone.montecarlo_radial1d(radial1d_mdl)
    spectrum = synspec.get_lambda_spec(out_nu, out_energy, 500*1e-8, 20000*1e-8, samples=1000)
    return out_nu, out_energy, jsestimator, nubarestimator, line_id_in, line_id_out, radial1d_mdl

#checking if kurucz_atom.h5 exists
kurucz_h5 = os.path.expanduser('~/.tardis/kurucz_atom.h5')
if not os.path.exists(kurucz_h5):
    raise IOError('kurucz_atom.h5 not found in %s' % kurucz_h5)


fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

config_object = config_reader.read_config('tardis_config_scatter.ini')

config_object.plasma_type='nebular'
config_object.line_interaction_type='macroatom'
config_object.single_run_packets = int(1e5)

atom_data = atomic.AtomData.from_hdf5(kurucz_h5, use_macro_atom=True, use_zeta_data=True)

out_nu, out_energy, jsestimator, nubarestimator, line_id_in, line_id_out, radial1d_mdl = plot_tardis()

emission_filter = out_nu >= 0
input = synspec.get_lambda_spec(radial1d_mdl.packet_src.packet_nus, radial1d_mdl.packet_src.packet_energies, 500*1e-8, 20000*1e-8, samples=1000)

line_id_in = np.array(line_id_in, dtype=int)
line_id_out = np.array(line_id_out, dtype=int)


line_in_atom = radial1d_mdl.atom_data.lines['atomic_number'].values[line_id_in]
line_out_atom = radial1d_mdl.atom_data.lines['atomic_number'].values[line_id_out]
assert all(line_in_atom == line_out_atom)
print "checked atoms"

line_in_ion = radial1d_mdl.atom_data.lines['ion_number'].values[line_id_in]
line_out_ion = radial1d_mdl.atom_data.lines['ion_number'].values[line_id_out]
assert all(line_in_ion == line_out_ion)
print "checked ions"


line_in_level_lower = radial1d_mdl.atom_data.lines['level_id_lower'].values[line_id_in]
line_out_level_upper = radial1d_mdl.atom_data.lines['level_id_upper'].values[line_id_out]








