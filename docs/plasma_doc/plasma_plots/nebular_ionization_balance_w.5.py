import os

from pylab import *
from tardis import atomic, plasma

atom_fname = os.path.expanduser('~/.tardis/si_kurucz.h5')

atom_data = atomic.AtomData.from_hdf5(atom_fname)
atom_data.prepare_atom_data([14], 'scatter')

nebular_plasma = plasma.NebularPlasma.from_abundance(10000, 0.5, {'Si': 1}, 1e-13, atom_data, 10.)
siI = []
siII = []
siIII = []
siIV = []
t_rads = linspace(5000, 20000, 100)
for t_rad in t_rads:
    nebular_plasma.update_radiationfield(t_rad, w=1.0)
    si_number_density = nebular_plasma.number_density.get_value(14)
    ion_density = nebular_plasma.ion_populations / si_number_density
    siI.append(ion_density.get_value((14, 0)))
    siII.append(ion_density.get_value((14, 1)))
    siIII.append(ion_density.get_value((14, 2)))
    siIV.append(ion_density.get_value((14, 3)))

plot(t_rads, siI, 'b-', label='Si I W=1.0')
plot(t_rads, siII, 'g-', label='Si II W=1.0')
plot(t_rads, siIII, 'r-', label='Si III W=1.0')
plot(t_rads, siIV, 'k-', label='Si IV W=1.0')

siI = []
siII = []
siIII = []
siIV = []
t_rads = linspace(5000, 20000, 100)
for t_rad in t_rads:
    nebular_plasma.update_radiationfield(t_rad, w=0.5)
    si_number_density = nebular_plasma.number_density.get_value(14)
    ion_density = nebular_plasma.ion_populations / si_number_density
    siI.append(ion_density.get_value((14, 0)))
    siII.append(ion_density.get_value((14, 1)))
    siIII.append(ion_density.get_value((14, 2)))
    siIV.append(ion_density.get_value((14, 3)))

plot(t_rads, siI, 'b--', label='Si I W=0.5')
plot(t_rads, siII, 'g--', label='Si II W=0.5')
plot(t_rads, siIII, 'r--', label='Si III W=0.5')
plot(t_rads, siIV, 'k--', label='Si IV W=0.5')

xlabel('T [K]')
ylabel('Number Density Fraction')

legend()
show()