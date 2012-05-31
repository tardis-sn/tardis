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

Bnu = lambda nu, t: (2 * constants.h * nu**3 / constants.c**2) / (np.exp(((constants.h * nu) / (constants.kb * t))) - 1)

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
    
def read_line_list(conn, atoms=None, symbol2z=None, max_atom=50, max_ion=50):
    raw_data = []
    
    select_stmt = """SELECT
                        lines.id, wl, loggf, g_lower, g_upper, e_lower,
                        e_upper, level_id_lower, level_id_upper,
                        global_level_id_lower, global_level_id_upper,
                        lines.atom, lines.ion, metastable
                    FROM
                        lines
                    inner join
                        levels
                    ON
                        global_level_id_lower=levels.id
                    
                    where
                        lines.atom < %(max_atom)d
                    and
                        lines.ion <  %(max_ion)d
                    %(select_atom_stmt)s
                    order by lines.id
                """ 
    
    if atoms is not None:
        if symbol2z is None:
            symbol2z = initialize.read_symbol2z()
            atom_list = ','.join([str(symbol2z[atom]) for atom in atoms])
            select_atom_stmt = """ and atom in (%s)""" % atom_list
    else:
        select_atom_stmt = ''
        
    line_select_stmt = select_stmt % dict(max_atom=max_atom, max_ion=max_ion, select_atom_stmt=select_atom_stmt)
    
    if sqlparse_available:
        print sqlparse.format(line_select_stmt, reindent=True)
    else:
        print line_select_stmt

    
    curs = conn.execute(line_select_stmt)

    
    for id, wl, loggf, g_lower, g_upper, e_lower, e_upper, level_id_lower, level_id_upper, global_level_id_lower, global_level_id_upper, atom, ion, metastable in curs:
        gf = 10**loggf
        f_lu = gf / g_lower
        f_ul = gf / g_upper
        raw_data.append((id, wl, g_lower, g_upper, f_lu, f_ul, e_lower, e_upper, level_id_lower, level_id_upper, global_level_id_lower, global_level_id_upper, atom, ion, metastable))
    
    line_list = np.array(raw_data, dtype = [('id', np.int64), ('wl', np.float64),
        ('g_lower', np.int64), ('g_upper', np.int64),
        ('f_lu', np.float64), ('f_ul', np.float64),
        ('e_lower', np.float64), ('e_upper', np.float64),
        ('level_id_lower', np.int64), ('level_id_upper', np.int64),
        ('global_level_id_lower', np.int64), ('global_level_id_upper', np.int64),
        ('atom', np.int64), ('ion', np.int64), ('metastable', np.bool)])
    
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
        self.target_level_id_down = np.array(reference_down)
        self.target_line_id_down = np.array(line_id_down)
        self.p_internal_down = np.array(p_internal_down)
        self.p_emission_down = np.array(p_emission_down)
        self.count_up = np.array(count_up)
        self.target_level_id_up = np.array(reference_up)
        self.target_line_id_up = np.array(line_id_up)
        self.p_internal_up = np.array(p_internal_up)
        
        self.count_total = 2 * self.count_down + self.count_up
        self.count_total_sum = np.sum(self.count_total)
        self.level_references = np.hstack(([0], np.cumsum(self.count_total)[:-1]))
    
    def merge_arrays(self):
        self.p_total = -1 * np.ones((self.count_total_sum), dtype=np.float64)
        self.transition_type_total = -1 * np.ones((self.count_total_sum), dtype=np.int64)
        self.target_line_total = -1 * np.ones((self.count_total_sum), dtype=np.int64)
        self.target_level_total = -1 * np.ones((self.count_total_sum), dtype=np.int64)
        for i, (ref, c_down, c_up, c_total) in enumerate(zip(self.level_references, self.count_down, self.count_up, self.count_total)):
            p_level_total = np.hstack((self.p_emission_down[i], self.p_internal_down[i], self.p_internal_up[i]))
            sort_probabilities = np.argsort(p_level_total)[::-1]
            transition_type_total = np.hstack((np.ones(c_down), np.zeros(c_down + c_up)))
            target_line_total = np.hstack((self.target_line_id_down[i], self.target_line_id_down[i], self.target_line_id_up[i]))
            target_level_total = np.hstack((self.target_level_id_down[i], self.target_level_id_down[i], self.target_level_id_up[i]))
            
            self.p_total[ref:ref + c_total] = p_level_total[sort_probabilities]
            self.transition_type_total[ref:ref + c_total] = transition_type_total[sort_probabilities]
            self.target_line_total[ref:ref + c_total] = target_line_total[sort_probabilities]
            self.target_level_total[ref:ref + c_total] = target_level_total[sort_probabilities]
        
        
        self.target_line_total -= 1
        self.target_level_total -= 1
        
        
    def calculate_beta_sobolev(self, tau_sobolev):
        
        def safe_beta_sobolev(tau_sobolev):
            if tau_sobolev > 1e3:
                return 1 / tau_sobolev
            elif tau_sobolev < 1e-4:
                return 1 - 0.5 * tau_sobolev
            else:
                return (1 - np.exp(-tau_sobolev)) / tau_sobolev
        vec_safe_beta_sobolev = np.vectorize(safe_beta_sobolev)
        beta_sobolev = vec_safe_beta_sobolev(tau_sobolev)
        
        return beta_sobolev
    
    def calculate_jbar(self, tau_sobolev, nu, t_rad, w, beta_sobolev=None):
        if beta_sobolev is None: beta_sobolev = self.calculate_beta_sobolev(tau_sobolev)
        jbar = beta_sobolev * w * Bnu(nu, t_rad)
        return jbar
    
    def calculate_p_downs(self, beta_sobolev):
        
        def beta_update(p_internal_down, p_emission_down, target_line_id, target_level_id):
            p_down = np.hstack((p_internal_down * beta_sobolev[target_line_id - 1], p_emission_down * beta_sobolev[target_line_id - 1]))
            #type is emitting or not
            type_down = np.hstack((np.zeros_like(p_internal_down, dtype=np.int64),
                                   np.ones_like(p_emission_down, dtype=np.int64)))
            target_level_id_down = np.hstack((target_level_id, target_level_id))
            target_line_id_down = np.hstack((target_line_id, target_line_id))
            return p_down, type_down, target_level_id_down, target_line_id_down
        vec_beta_update = np.vectorize(beta_update, otypes = (np.ndarray, np.ndarray, np.ndarray, np.ndarray))
        return vec_beta_update(self.p_internal_down, self.p_emission_down, self.target_line_id_down, self.target_level_id_down)
    
    def calculate_p_ups(self, jbar):
        
        def jbar_update(p_internal_up, target_line_id):
            p_up = p_internal_up * jbar[target_line_id - 1]
            #type is emitting or not
            type_up = np.zeros_like(p_up, dtype=np.int64)
            return p_up, type_up
        
        vec_jbar_update = np.vectorize(jbar_update, otypes = (np.ndarray,np.ndarray))
        return vec_jbar_update(self.p_internal_up, self.target_line_id_up)
    
    def calculate_transition_probabilities(self, beta_sobolev, jbar):
        p_down, type_down, target_level_id_down, target_line_id_down = self.calculate_p_downs(beta_sobolev)
        p_up, type_up = self.calculate_p_ups(jbar)
        
        #merging probabilites
        def merge_probability (p_down_array, p_up_array,
                               type_down_array, type_up_array,
                               target_level_id_down_array, target_level_id_up_array,
                               target_line_id_down_array, target_line_id_up_array):
            p_array = np.hstack((p_up_array, p_down_array))
            norm_p_array = np.sum(p_array)
            if norm_p_array > 0.0:
                p_array = p_array / norm_p_array 
            return (p_array,
                    np.hstack((type_up_array, type_down_array)),
                    np.hstack((target_level_id_up_array, target_level_id_down_array)),
                    np.hstack((target_line_id_up_array, target_line_id_down_array)))
            
        vec_merge_probability = np.vectorize(merge_probability, otypes=(np.ndarray, np.ndarray, np.ndarray, np.ndarray))
        p_transition, type_transition, target_level_id, target_line_id = \
                                    vec_merge_probability(p_down, p_up,
                                                          type_down, type_up,
                                                          target_level_id_down, self.target_level_id_up,
                                                          target_line_id_down, self.target_line_id_up)
        p_transition_merged = np.hstack(p_transition)
        type_transition_merged = np.hstack(type_transition)
        target_level_id_merged = np.hstack(target_level_id) - 1
        target_line_id_merged = np.hstack(target_line_id) - 1
        no_probabilities = self.count_down * 2 + self.count_up
        unroll_reference = np.hstack(([0], np.cumsum(no_probabilities)[:-1]))
        return p_transition_merged, type_transition_merged, target_line_id_merged, target_level_id_merged, no_probabilities, unroll_reference
    
    