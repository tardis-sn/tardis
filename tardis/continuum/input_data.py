from astropy import constants as const
from astropy import units as units
import numpy as np


class ContinuumInputData(object):
    c_einstein = (4. * (np.pi * const.e.esu) ** 2 / (const.c.cgs * const.m_e.cgs)).value
    C0_ff = 1.426e-27  # in cgs units (see Osterbrock 1974)
    c_0_regemorter = 5.465e-11
    I_H = 13.598433770784 * units.eV.to(units.erg)

    def __init__(self, atom_data, plasma_array, ws):
        # Plasma quantities
        self.electron_densities = plasma_array.electron_densities.values
        self.t_electrons = plasma_array.t_electrons
        self.t_rads = plasma_array.t_rad
        self.link_t_rad_t_electron = plasma_array.link_t_rad_t_electron
        self.ion_number_density = plasma_array.ion_number_density
        # TODO: Replace level population with LTE level population
        self.lte_level_pop = plasma_array.level_number_density
        self.level_pop = plasma_array.level_number_density

        # Radiation field
        self.ws = ws

        # Atom data
        self.lines = atom_data.lines
        self.levels = atom_data.levels
        self.photoionization_data = atom_data.continuum_data.photoionization_data

        #
        self.macro_atom_references = atom_data.macro_atom_references
        self.macro_atom_data = atom_data.macro_atom_data
        self.macro_atom_continuum_data = atom_data.continuum_data.macro_atom_data

        #
        self.continuum_references = atom_data.continuum_data.continuum_references

        ##
        self.ion_charges = self._get_ion_charges()

        # Computed quantities
        self.nu_i = self._get_nu_i()


    def _get_nu_i(self):
        return self.photoionization_data.groupby(level=[0, 1, 2]).first().nu.values

    def _get_ion_charges(self):
        return self.ion_number_density.index.get_level_values(1).values