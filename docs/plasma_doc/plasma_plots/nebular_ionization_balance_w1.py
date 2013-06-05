import os

from pylab import *
from tardis import atomic, plasma

atom_fname = os.path.expanduser('~/.tardis/kurucz_atom.h5')

atom_data = atomic.AtomData.from_hdf5(atom_fname)
atom_data.prepare_atom_data([14], 'scatter')

nebular_plasma = plasma.NebularPlasma.from_abundance(1000, 1.0, {'Si': 1}, 1, atom_data, 10.)
siI = []
siII = []
siIII = []
siIV = []
t_rads = linspace(5000, 20000, 100)
for t_rad in t_rads:
    nebular_plasma.update_radiationfield(t_rad, w=1.)
    si_number_density = nebular_plasma.number_density.get_value(14)
    ion_density = nebular_plasma.ion_populations / si_number_density
    siI.append(ion_density.get_value((14, 0)))
    siII.append(ion_density.get_value((14, 1)))
    siIII.append(ion_density.get_value((14, 2)))
    siIV.append(ion_density.get_value((14, 3)))

plot(t_rads, siI, label='Si I')
plot(t_rads, siII, label='Si II')
plot(t_rads, siIII, label='Si III')
plot(t_rads, siIV, label='Si IV')

xlabel('T [K]')
ylabel('Number Density Fraction')
title('W=1.0')
legend()
show()