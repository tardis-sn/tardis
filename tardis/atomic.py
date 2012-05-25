# atomic model

from scipy import interpolate
import constants
import initialize
import numpy as np


class AtomModel(object):
    pass




class CompleteKuruczAtomModel(AtomModel):
    @classmethod
    def from_db(cls, conn, max_atom=30, max_ion=30):
        
        symbol2z = initialize.read_symbol2z()
        
        masses = initialize.read_atomic_data()['mass'][:max_atom]
        
        ionization_data = initialize.read_ionization_data()
        
        
        #reason for max_ion - 1: in energy level data there's unionized, once-ionized, twice-ionized, ...
        #in ionization_energies, there's only once_ionized, twice_ionized
        ionization_energy = np.zeros((max_ion - 1, max_atom))
        
        for atom, ion, ion_energy in ionization_data:
            if atom > max_atom or ion >= max_ion:
                continue
            ionization_energy[ion-1, atom-1] = ion_energy
        
        
        levels_energy, levels_g = read_kurucz_level_data_fromdb(conn, max_atom, max_ion)
        
        #factor zeta ML 1993
        recombination_coefficient_interp = read_recombination_coefficients_fromdb(conn, max_atom, max_ion)
        
        
        
        
        return cls(masses=masses, ionization_energy=ionization_energy,
                   levels_energy=levels_energy, levels_g=levels_g,
                   recombination_coefficient_interp=recombination_coefficient_interp,
                   max_atom=max_atom, max_ion=max_ion)
    
    def __init__(self, masses=None, ionization_energy=None,
                levels_energy=None,
                levels_g=None,
                recombination_coefficient_interp=None,
                max_atom=None, max_ion=None):
        self.masses = masses
        self.ionization_energy = ionization_energy
        self.levels_energy = levels_energy
        self.levels_g = levels_g
        #unclear switch variable names around?
        self.interpolate_recombination_coefficient = recombination_coefficient_interp
        self.max_atom = max_atom
        self.max_ion = max_ion
    
    def calculate_radfield_correction_factor(self, t_rad, t_electron, w, departure_coefficient=None, xi_threshold_species = (1, 19)):
    #factor delta ML 1993
        if departure_coefficient is None:
            departure_coefficient = 1 / float(w)
        delta = np.ones((self.max_ion - 1, self.max_atom))
        xi_threshold = self.ionization_energy[xi_threshold_species]
        
        #Formula 15 ML 1993
        threshold_filter = (self.ionization_energy <= xi_threshold) & (self.ionization_energy > 0)
        delta[threshold_filter] = (t_electron / (departure_coefficient * w * t_rad)) * \
            np.exp((delta[threshold_filter] / (constants.kbinev * t_rad)) \
            - (delta[threshold_filter] / (constants.kbinev * t_electron)))
        
        
        threshold_filter = (self.ionization_energy > xi_threshold) & (self.ionization_energy > 0)
        #Formula 20 ML 1993
        delta[self.ionization_energy > xi_threshold] = 1 - \
            np.exp((delta[threshold_filter] / (constants.kbinev * t_rad)) \
            - (xi_threshold / (constants.kbinev * t_rad))) + \
            (t_electron / (departure_coefficient * w * t_rad)) * \
            np.exp((delta[threshold_filter] / (constants.kbinev * t_rad)) \
            - (delta[threshold_filter] / (constants.kbinev * t_electron)))
        
        return delta
        

def read_recombination_coefficients_fromdb(conn, max_atom=30, max_ion=None):
    
    t_rads = np.linspace(2000, 50000, 10)
    recombination_coefficients = np.ones((max_ion - 1, max_atom, 10)) * 0.5
    interpolator = interpolate.interp1d(t_rads, recombination_coefficients, kind='linear')
    return interpolator
    
def read_kurucz_level_data_fromdb(conn, max_atom=30, max_ion=None):
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
