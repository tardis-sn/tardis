#calculation for line interaction
import constants
import numpy as np
import initialize

try:
    import sqlparse
    sqlparse_available = True
except ImportError:
    sqlparse_available = False
    
def get_tau_line(nu, f_line, dens_line, t_exp):
    tau_line = constants.sobolev_coeff * f_line * (constants.c / nu) * dens_line * t_exp
    return tau_line

def get_r_sobolev(r, mu, d_line):
    return np.sqrt(r**2 + d_line**2 + 2 * r * d_line * mu)
    
    
def read_line_list(conn, atoms=None, symbol2z=None):
    
            
    
    
    raw_data = []
    
    select_stmt = """SELECT
                        wl, loggf, g_lower, g_upper, e_lower,
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
        curs = conn.execute(select_stmt % atom_list)
        print atom_list
    else:
        curs = conn.execute(select_stmt)
    
    
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