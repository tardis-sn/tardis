#Calculations of the Plasma conditions

import sqlite3
import constants
import math
import numpy as np
import scipy.sparse
import pdb
kbinev = constants.kb * constants.erg2ev
h = 4.135667516e-15 #eV * s
h_cgs = 6.62606957e-27
u = 1.66053886e-24 # atomic mass unit in g
me = 9.10938188e-28 #grams



class Plasma(object):
    pass

class LTEPlasma(Plasma):
    
    

def get_atomic_data(conn):
    data = conn.execute('select atom, symbol, weight from atoms').fetchall()
    return np.array(data,
        dtype=[('atom', np.int64), ('symbol', 'S2'), ('weight', np.float64)])
    
def make_levels_table(conn):
    conn.execute('create temporary table levels as select elem, ion, e_upper * ? * ? as energy, j_upper as j, label_upper as label from lines union select elem, ion, e_lower, j_lower, label_lower from lines', (constants.c, h))
    conn.execute('create index idx_all on levels(elem, ion, energy, j, label)')

def get_ionize_data(conn, max_atom=30, max_ion=None):
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
    

def get_level_data(conn, max_atom=30, max_ion=None):
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
                    ion <= ?
            order by
                atom, ion, energy"""
                    
    curs = conn.execute(level_select_stmt, (max_ion, max_atom))
    energy_data = np.zeros((max_ion + 1, max_atom), dtype='object')
    g_data = np.zeros((max_ion + 1, max_atom), dtype='object')
    
    
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
    
def calculate_partition_functions(energy_data, g_data, beta):
    z_func = lambda energy, g: np.sum(g*np.exp(-beta*energy))
    vector_z_func = np.vectorize(z_func, otypes=[np.float64])
    return vector_z_func(energy_data, g_data)


def calculate_phis(partition_functions, ionize_data, beta):
    #calculating ge = 2/(Lambda^3)
    ge = 2 / (np.sqrt(h_cgs**2 / (2 * np.pi * me * (1 / (beta*constants.erg2ev)))))**3
    
    partition_fractions = ge * partition_functions[1:]/partition_functions[:-1]
    partition_fractions[np.isnan(partition_fractions)] = 0.0
    partition_fractions[np.isinf(partition_fractions)] = 0.0
    
    phis = partition_fractions * np.exp(-beta * ionize_data)
    
    return phis
    
    
def calculate_single_ion_populations(phis, atom_number_density, electron_density=None, max_atom=99):
    #N1 is ground state
    #N(fe) = N1 + N2 + .. = N1 + (N2/N1)*N1 + (N3/N2)*(N2/N1)*N1 + ... = N1(1+ N2/N1+...)
    if electron_density == None:
        electron_density = np.sum(atom_number_density)
    
    ion_fraction_prod = np.cumprod(phis / electron_density, axis = 0) # (N2/N1, N3/N2 * N2/N1, ...)
    ion_fraction_sum = 1 + np.sum(ion_fraction_prod, axis = 0)
    N1 = atom_number_density / ion_fraction_sum
    #Further Ns
    Nn = N1 * ion_fraction_prod
    new_electron_density = np.sum(Nn * (np.arange(max_atom)+1).reshape((max_atom,1)))
    return np.vstack((N1, Nn)), new_electron_density

def calculate_ion_populations(beta, W, atom_number_density,
                            partition_functions, atomic_data, ionize_data,
                            max_atom=30, max_ion=30):
    #partition_functions = calculate_partition_functions(energy_data, g_data, beta)
    phis = calculate_phis(partition_functions, ionize_data, beta)
    
    #first estimate
    old_electron_density = np.sum(atom_number_density)
    while True:
        ion_density, electron_density = \
            calculate_single_ion_populations(phis, atom_number_density, old_electron_density, max_atom=max_atom)
        
        if abs(electron_density / old_electron_density - 1) < 0.05: break
        #print "ion population calculation %s" % (abs(electron_density / old_electron_density - 1),)
        old_electron_density = 0.5*(old_electron_density + electron_density)
    return ion_density, electron_density
        
def calculate_atom_number_density(named_abundances, density, atomic_data, z2symbol, max_atom=30):
    number_density = np.zeros(max_atom)
    assert sum(named_abundances.values()) == 1.
    
    for atom in xrange(max_atom):
        symbol = z2symbol[atom + 1]
        abundance = named_abundances.get(symbol, 0.)
        number_density[atom] = (density * abundance) / (atomic_data['weight'][atom - 1] * u)
    return number_density


