#calculation for line interaction
import constants
import numpy as np
import initialize
import sqlite3

try:
    import sqlparse
    sqlparse_available = True
except ImportError:
    sqlparse_available = False
    
def get_tau_line(nu, f_line, dens_line, t_exp):
    tau_line = constants.sobolev_coeff * f_line * (constants.c / nu) * dens_line * t_exp
    return tau_line

def convert_int_ndarray(sqlite_binary):
    if sqlite_binary == '-1':
        return np.array([], dtype=np.int64)
    else:
        return np.frombuffer(sqlite_binary, dtype=np.int64)

def convert_float_ndarray(sqlite_binary):
    if sqlite_binary == '-1.0':
        return np.array([], dtype=np.float64)
    else:
        return np.frombuffer(sqlite_binary, dtype=np.float64)

sqlite3.register_converter('int_ndarray', convert_int_ndarray)
sqlite3.register_converter('float_ndarray', convert_float_ndarray)
    
def read_line_list(conn, atoms=None, symbol2z=None):
    
            
    
    
    raw_data = []
    
    select_stmt = """SELECT
                        id, wl, loggf, g_lower, g_upper, e_lower,
                        e_upper, level_id_lower, level_id_upper,
                        atom, ion
                    FROM
                        lines
                    
                """
    
    if atoms is not None:
        if symbol2z is None:
            symbol2z = initialize.read_symbol2z()
            select_stmt += """ where atom in (%s)"""
            
        atom_list = ','.join([str(symbol2z[atom]) for atom in atoms])
        curs = conn.execute(select_stmt % atom_list + 'order by wl')
        print atom_list
    else:
        curs = conn.execute(select_stmt)
    
    
    for id, wl, loggf, g_lower, g_upper, e_lower, e_upper, level_id_lower, level_id_upper, atom, ion in curs:
        gf = 10**loggf
        f_lu = gf / g_lower
        f_ul = gf / g_upper
        raw_data.append((id, wl, g_lower, g_upper, f_lu, f_ul, e_lower, e_upper, level_id_lower, level_id_upper, atom, ion))
    
    line_list = np.array(raw_data, dtype = [('id', np.int64), ('wl', np.float64),
        ('g_lower', np.int64), ('g_upper', np.int64),
        ('f_lu', np.float64), ('f_ul', np.float64),
        ('e_lower', np.float64), ('e_upper', np.float64),
        ('level_id_lower', np.int64), ('level_id_upper', np.int64),
        ('atom', np.int64), ('ion', np.int64)])
    
    return line_list



class MacroAtomData(object):
    pass

class SimpleMacroAtomData(MacroAtomData):
    
    @classmethod
    def fromdb(cls, conn):
        curs = conn.execute("""
        SELECT
            count_down,
            reference_down,
            line_id_down,
            p_internal_down,
            p_emission_down,
            count_up,
            reference_up,
            line_id_up,
            p_internal_up
        FROM
            macro_atom
        ORDER BY
            id""")
        
        (count_down, reference_down, line_id_down, p_internal_down, p_emission_down,
         count_up, reference_up, line_id_up, p_internal_up) = zip(*curs.fetchall())
        
        return cls(count_down, reference_down, line_id_down, p_internal_down, p_emission_down,
         count_up, reference_up, line_id_up, p_internal_up)
        
        
    def __init__(self,
                count_down, reference_down, line_id_down,
                p_internal_down, p_emission_down,
                count_up, reference_up, line_id_up, p_internal_up):
        self.count_down = np.array(count_down)
        self.reference_down = np.array(reference_down)
        self.line_id_down = np.array(line_id_down)
        self.p_internal_down = np.array(p_internal_down)
        self.p_emission_down = np.array(p_emission_down)
        self.count_up = np.array(count_up)
        self.reference_up = np.array(reference_up)
        self.line_id_up = np.array(line_id_up)
        self.p_internal_up = np.array(p_internal_up)
        
    
def compile_sobolev_tau(line_list, partition_functions, ion_population, t_exp):
    tau_sobolev = []
    wl = line_list['wl'] * 1e-7
    C = 1 #supposed to be (pi*e**2)/(m_e * c)
    Z = partition_functions[line_list['ion'], line_list['atom']-1]
    g_lower = line_list['g_lower']
    e_lower = line_list['e_lower']
    n_lower = (g_lower / Z) * np.exp(-beta * e_lower) * ion_populations[line_list['ion'], line_list['atom']-1]
    tau_sobolev = C * line_list['f_lu'] * wl * t_exp * n_lower
    return tau_sobolev