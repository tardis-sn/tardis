#Calculations of the Plasma conditions

import sqlite3
import constants
import math
import numpy as np

kbinev = constants.kb * constants.erg2ev
h = 4.135667516e-15 #eV * s
class PartitionLTE:
    def __init__(self):
        self.Z = 0
    
    def step(self, energy, j, beta):
        energy *= constants.c * h
        self.Z += (2*j + 1) * math.exp(-beta * energy)

    def finalize(self):
        return self.Z


class TestAggregate:
    def __init__(self):
        self.Z = 0
    
    def step(self, energy, j, beta):
        print energy
        self.Z += (2*j+1)*math.exp(-beta * energy)
    
    def finalize(self):
        return self.Z


def make_levels_table(conn):
    conn.execute('create temporary table levels as select elem, ion, e_upper * ? * ? as energy, j_upper as j, label_upper as label from lines union select elem, ion, e_lower, j_lower, label_lower from lines', (constants.c, h))
    conn.execute('create index idx_all on levels(elem, ion, energy, j, label)')

def calculate_partition(T, conn):
    beta = 1 / (T * constants.kb * constants.erg2ev)
    #conn.create_aggregate('lte_partition', 3, PartitionLTE)
    partition_function_sql = """select element, ion_number, ionize_ev, partition
                                from
                                    (select levels.elem as element, levels.ion as ion_number, lte_partition(energy, j, ?) as partition
                                    from levels
                                    group by levels.elem, levels.ion)
                                left join ionization
                                    on ionization.atom=element and ionization.ion = ion_number
                                order by element, ion_number
                                """
    #partition_function_sql = """select levels.elem as element, levels.ion, lte_partition(energy, j, ?) as partition
    #                                from levels
    #                                group by levels.elem, levels.ion    
    #                            """
    partition_functions = np.array(conn.execute(partition_function_sql, (beta,)).fetchall())
    return partition_functions

def calculate_phis(atom, ion_upper, beta, partition_functions):
    upper_idx = np.where(np.logical_and(partition_functions[:,0] == atom, partition_functions[:,1] == ion))[0]
    lower_idx = np.where(np.logical_and(partition_functions[:,0] == atom, partition_functions[:,1] == ion-1))[0]
    
    #comment out for performance maybe
    assert upper_idx.size == 1 #asserting that there's only one entry for the upper ion level
    assert lower_idx.size == 1 #asserting that there's only one entry for the lower ion level
    
    
    
    z_upper, z_lower = partition_functions[upper_idx][0][-1], partition_functions[lower_idx][0][-1]
    ionization_energy = partition_functions[upper_idx][0][-2]
    
    phi = (z_upper/z_lower) * math.exp(-beta*ionization_energy)
    return phi
    
    
    