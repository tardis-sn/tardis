#Calculations of the Plasma conditions

import sqlite3
import constants
import math
import numpy as np
import scipy.sparse
import pdb
import initialize

kbinev = constants.kb * constants.erg2ev
h = 4.135667516e-15 #eV * s
h_cgs = 6.62606957e-27
u = 1.66053886e-24 # atomic mass unit in g
me = 9.10938188e-28 #grams

    

class Plasma(object):
    pass

class LTEPlasma(Plasma):
    
    @classmethod
    def from_db(cls, named_abundances, density, conn, max_atom=30, max_ion=30):
        
        
        
        symbol2z = initialize.read_symbol2z()
        
        #converting the abundances dictionary dict to an array
        abundances = np.zeros(max_atom)
        
        for symbol in named_abundances:
            abundances[symbol2z[symbol]-1] = named_abundances[symbol]
        
        
        masses = initialize.read_atomic_data()['mass'][:max_atom]
        
        ionization_data = initialize.read_ionization_data()
        
        
        #reason for max_ion - 1: in energy level data there's unionized, once-ionized, twice-ionized, ...
        #in ionization_energies, there's only once_ionized, twice_ionized
        ionize_energy = np.zeros((max_ion - 1, max_atom))
        
        for atom, ion, ion_energy in ionization_data:
            if atom > max_atom or ion >= max_ion:
                continue
            ionize_energy[ion-1, atom-1] = ion_energy
        
        
        levels_energy, levels_g = read_level_data(conn, max_atom, max_ion)
        
        return cls(abundances, density,
                   masses=masses, ionize_energy=ionize_energy,
                   levels_energy=levels_energy, levels_g=levels_g,
                   max_atom=max_atom, max_ion=max_ion)
        
            
    def __init__(self, abundances, density,
                masses=None, ionize_energy=None,
                levels_energy=None,
                levels_g=None, 
                max_atom=None, max_ion=None):
        """
        ionize_energy
        -------------
        ndarray with row being ion, column atom
        
        """
        self.masses = masses
        self.ionize_energy = ionize_energy
        self.levels_energy = levels_energy
        self.levels_g = levels_g
        self.max_atom = max_atom
        self.max_ion = max_ion
        
        self.atom_number_density = self._calculate_atom_number_density(abundances, density)
        self.electron_density = np.sum(self.atom_number_density)


    def update_radiationfield(self, t_rad):
        self.t_rad = t_rad
        self.beta = 1 / (t_rad * constants.kbinev)
        self.partition_functions = self._calculate_partition_functions()
        self.ion_number_density, self.electron_density = self._calculate_ion_populations()
        
    
    
    def calculate_tau_sobolev(self, line_list, time_exp):
        wl = line_list['wl'] * 1e-7
        C = (np.pi * constants.e**2) / (constants.me * constants.c) #supposed to be (pi*e**2)/(m_e * c)
        Z = self.partition_functions[line_list['ion'], line_list['atom']-1]
        g_lower = line_list['g_lower']
        e_lower = line_list['e_lower']
        n_lower = (g_lower / Z) * np.exp(-self.beta * e_lower) * self.ion_number_density[line_list['ion'], line_list['atom']-1]
        tau_sobolev = C * line_list['f_lu'] * wl * time_exp * n_lower
        return tau_sobolev
    
    def _calculate_atom_number_density(self, abundances, density):
        #TODO float comparison problematic
        assert sum(abundances) == 1.
        number_density = (density * abundances) / (self.masses * u)
    
        return number_density
    
    def _calculate_partition_functions(self):
        z_func = lambda energy, g: np.sum(g*np.exp(-self.beta*energy))
        vector_z_func = np.vectorize(z_func, otypes=[np.float64])
        return vector_z_func(self.levels_energy, self.levels_g)

    def _calculate_phis(self):
        #calculating ge = 2/(Lambda^3)
        ge = 2 / (np.sqrt(h_cgs**2 / (2 * np.pi * me * (1 / (self.beta*constants.erg2ev)))))**3
        
        partition_fractions = ge * self.partition_functions[1:] / self.partition_functions[:-1]
        partition_fractions[np.isnan(partition_fractions)] = 0.0
        partition_fractions[np.isinf(partition_fractions)] = 0.0
        #phi = (n_j+1 * ne / nj)
        phis = partition_fractions * np.exp(-self.beta * self.ionize_energy)
        return phis

    def _calculate_single_ion_populations(self, phis):
        #N1 is ground state
        #N(fe) = N1 + N2 + .. = N1 + (N2/N1)*N1 + (N3/N2)*(N2/N1)*N1 + ... = N1(1+ N2/N1+...)
        
        ion_fraction_prod = np.cumprod(phis / self.electron_density, axis = 0) # (N2/N1, N3/N2 * N2/N1, ...)
        ion_fraction_sum = 1 + np.sum(ion_fraction_prod, axis = 0)
        N1 = self.atom_number_density / ion_fraction_sum
        #Further Ns
        Nn = N1 * ion_fraction_prod
        new_electron_density = np.sum(Nn * (np.arange(1,self.max_ion).reshape((self.max_ion - 1,1))))
        return np.vstack((N1, Nn)), new_electron_density

    def _calculate_ion_populations(self):
        #partition_functions = calculate_partition_functions(energy_data, g_data, beta)
        phis = self._calculate_phis()
    
        #first estimate
        electron_density = np.sum(self.atom_number_density)
        old_electron_density = self.electron_density
        while True:
            ion_density, electron_density = \
                self._calculate_single_ion_populations(phis)
            
            if abs(electron_density / old_electron_density - 1) < 0.05: break
            
            old_electron_density = 0.5*(old_electron_density + electron_density)
        return ion_density, electron_density
    
    
    
    
class NebularPlasma(LTEPlasma):
    
    @classmethod
    def from_db(cls, named_abundances, density, conn, max_atom=30, max_ion=30):
        symbol2z = initialize.read_symbol2z()
        
        #converting the abundances dictionary dict to an array
        abundances = np.zeros(max_atom)
        
        for symbol in named_abundances:
            abundances[symbol2z[symbol]-1] = named_abundances[symbol]
        
        
        masses = initialize.read_atomic_data()['mass'][:max_atom]
        
        ionization_data = initialize.read_ionization_data()
        
        
        #reason for max_ion - 1: in energy level data there's unionized, once-ionized, twice-ionized, ...
        #in ionization_energies, there's only once_ionized, twice_ionized
        ionize_energy = np.zeros((max_ion - 1, max_atom))
        
        for atom, ion, ion_energy in ionization_data:
            if atom > max_atom or ion >= max_ion:
                continue
            ionize_energy[ion-1, atom-1] = ion_energy
        
        
        levels_energy, levels_g = read_level_data(conn, max_atom, max_ion)
        
        return cls(abundances, density,
                   masses=masses, ionize_energy=ionize_energy,
                   levels_energy=levels_energy, levels_g=levels_g,
                   max_atom=max_atom, max_ion=max_ion)
    
    def __init__(self, abundances, density,
                masses=None, ionize_energy=None,
                levels_energy=None,
                levels_g=None,
                max_atom=None, max_ion=None):

        LTEPlasma.__init__(self, abundances, density,
                masses=masses, ionize_energy=ionize_energy,
                levels_energy=levels_energy,
                levels_g=levels_g,
                max_atom=max_atom, max_ion=max_ion)
    
    def update_radiationfield(self, t_rad, w, t_electron=None):
        self.t_rad = t_rad
        self.w = w
        
        if t_electron is None:
            self.t_electron = 0.9 * self.t_rad
            
        self.beta = 1 / (t_rad * constants.kbinev)
        self.partition_functions = self._calculate_partition_functions()
        self.ion_number_density, self.electron_density = self._calculate_ion_populations()
        
        
    def _calculate_phis(self):
        #calculating ge = 2/(Lambda^3)
        ge = 2 / (np.sqrt(h_cgs**2 / (2 * np.pi * me * (1 / (self.beta*constants.erg2ev)))))**3
        
        partition_fractions = ge * self.partition_functions[1:] / self.partition_functions[:-1]
        partition_fractions[np.isnan(partition_fractions)] = 0.0
        partition_fractions[np.isinf(partition_fractions)] = 0.0
        #phi = (n_j+1 * ne / nj)
        phis = partition_fractions * np.exp(-self.beta * self.ionize_energy)
        phis = self.w * (self.t_electron/ self.t_rad)**.5 * phis
        return phis

    



def read_line_list(conn):
    
    raw_data = []
    curs = conn.execute('select wl, loggf, g_lower, g_upper, e_lower, e_upper, level_id_lower, level_id_upper, atom, ion from lines')
    for wl, loggf, g_lower, g_upper, e_lower, e_upper, level_id_lower, level_id_upper, atom, ion in curs:
        gf = 10**loggf
        f_lu = gf / g_lower
        f_ul = gf / g_upper
        raw_data.append((wl, g_lower, g_upper, f_lu, f_ul, e_lower, e_upper, level_id_lower, level_id_upper, atom, ion))
    
    line_list = np.array(raw_data, dtype = [('wl', np.float64),
        ('g_lower', np.int64), ('g_upper', np.int64),
        ('f_lu', np.float64), ('f_ul', np.float64),
        ('e_lower', np.float64), ('e_upper', np.float64),
        ('level_id_lower', np.int64), ('level_id_upper', np.int64),
        ('atom', np.int64), ('ion', np.int64)])
    
    return line_list


        
def make_levels_table(conn):
    conn.execute('create temporary table levels as select elem, ion, e_upper * ? * ? as energy, j_upper as j, label_upper as label from lines union select elem, ion, e_lower, j_lower, label_lower from lines', (constants.c, h))
    conn.execute('create index idx_all on levels(elem, ion, energy, j, label)')



    

def read_ionize_data_from_db(conn, max_atom=30, max_ion=None):
    if max_ion == None: max_ion = max_atom
    ionize_data = np.zeros((max_ion, max_atom))
    ionize_select_stmt = """select
                    atom, ion, ionize_ev
                from
                    ionization
                where
                    atom <= ?
                and
                    ion <= ?"""
    
    curs = conn.execute(ionize_select_stmt, (max_atom, max_ion))
    
    for atom, ion, ionize in curs:
        ionize_data[ion-1, atom-1] = ionize
    
    return ionize_data
    

def read_level_data(conn, max_atom=30, max_ion=None):
    #Constructing Matrix with atoms columns and ions rows
    #dtype is object and the cells will contain arrays with the energy levels
    if max_ion == None:
        max_ion = max_atom
    level_select_stmt = """select
                atom, ion, energy, g, level_id
            from
                levels
            where
                    atom <= ?
                and
                    ion < ?
            order by
                atom, ion, energy"""
                    
    curs = conn.execute(level_select_stmt, (max_ion, max_atom))
    energy_data = np.zeros((max_ion, max_atom), dtype='object')
    g_data = np.zeros((max_ion, max_atom), dtype='object')
    
    
    old_elem = None
    old_ion = None
    
    for elem, ion, energy, g, levelid in curs:
        if elem == old_elem and ion == old_ion:
            energy_data[ion, elem - 1] = np.append(energy_data[ion, elem - 1], energy)
            g_data[ion, elem - 1] = np.append(g_data[ion, elem - 1], g)
            
        else:
            old_elem = elem
            old_ion = ion
            energy_data[ion, elem - 1] = np.array([energy])
            g_data[ion, elem - 1] = np.array([g])
            
    return energy_data, g_data

def get_transition_data(conn, max_atom=99):
    curs = conn.execute('select atom, ion, loggf, level_id_upper, levelid_lower from lines')
    #continue writing for branching


