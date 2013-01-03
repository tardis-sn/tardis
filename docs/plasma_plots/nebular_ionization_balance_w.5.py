import os

from pylab import *
from tardis import atomic, plasma, model_radial_oned

atom_fname = os.path.expanduser('~/.tardis/kurucz_atom.h5')

atom_data = atomic.AtomData.from_hdf5(atom_fname, use_zeta_data=True)
number_density = model_radial_oned.calculate_atom_number_densities(atom_data, {'Si':1}, density=1e-13)
atom_data.prepare_atom_data(set(number_density.index.values.astype(int)), max_ion_number=3)
nebular_plasma = plasma.NebularPlasma(number_density, atom_data, max_ion_number=3)
siI = []
siII = []
siIII = []
siIV = []
t_rads = linspace(5000, 20000, 100)
for t_rad in t_rads:
    nebular_plasma.update_radiationfield(t_rad, w=.5)
    si_number_density = number_density.get_value(14, 'number_density')
    ion_density = nebular_plasma.ion_number_density / si_number_density
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
title('W=0.5')
legend()
show()