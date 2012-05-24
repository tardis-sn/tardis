#testing the code and new functions before they are sorted into a file

import initialize
import plasma
import line
import montecarlo
import c_montecarlo
import model
import photon
import constants
import synspec
import numpy as np
import cProfile

def run_simple_test(cfname, conn, no_of_packets=int(1e5)):
    T, t_exp, vph, dens, named_abundances = initialize.get_simple_config(cfname)
    beta = T * constants.kbinev
    
    r_inner = t_exp * vph
    r_outer = 1.2 * r_inner
    v_outer = vph * (r_outer/r_inner)
    
    W=1.
    symbol2z = initialize.get_symbol2z(conn)
    z2symbol = initialize.get_z2symbol(conn)
    atomic_data = plasma.get_atomic_data(conn)
    ionize_data = plasma.get_ionize_data(conn)
    e_data, g_data = plasma.get_level_data(conn)
    
    atom_number_densities = plasma.calculate_atom_number_density(named_abundances, dens, atomic_data, z2symbol, max_atom=30)
    partition_functions = plasma.calculate_partition_functions(e_data, g_data, beta)
    
    ion_number_densities, ne = plasma.calculate_ion_populations(beta,
                                                         W,
                                                         atom_number_densities,
                                                         partition_functions=partition_functions,
                                                         atomic_data=atomic_data,
                                                         ionize_data=ionize_data,
                                                         max_atom=30,
                                                         max_ion=30)
    line_list = line.read_line_list(conn)
    
    tau_sobolevs = line.compile_tau_sobolev(line_list, beta, partition_functions, ion_number_densities, t_exp)
    
    #testing without lines
    #tau_sobolevs = np.ones_like(tau_sobolevs) * 500.
    
    #return line_list, tau_sobolev, ion_number_densities, atom_number_densities, partition_functions
    simpleModel = model.OneZoneBaseModel(r_inner, r_outer, v_outer, ne)
    nus = photon.random_blackbody_nu(T, nu_range=(1e8*constants.c / 2000, 1e8*constants.c/ 9000), size=no_of_packets)
    mus = np.sqrt(np.random.random(no_of_packets))
    
    wl = line_list['wl']*1e-8
    line_list_nu = constants.c / wl
    #nu, energy = montecarlo.runSimpleOneD(photons, mus, line_list['wl']*1e-8, tau_sobolev, simpleModel)
    #cProfile.runctx("""nu, energy = c_montecarlo.run_simple_oned(photons, mus, line_list_nu, tau_sobolevs,
#                    r_inner, r_outer, vph, ne)""", globals(), locals(), "c_run_test.prof")
    #return None
    print r_inner, r_outer
    nu, energy, j_estimator, nubar_estimator = c_montecarlo.run_simple_oned(nus, mus, line_list_nu, tau_sobolevs,
                    r_inner, r_outer, vph, ne)
    nu_filter = nu > 0
    #energy = energy[nu_filter]
    #nu = nu[nu_filter]
    print "j_est %.2f nubar_est %.2f" % (j_estimator, nubar_estimator)
    print "ratio %.5f" % (j_estimator/nubar_estimator)
    return nu, energy, j_estimator, nubar_estimator
    
    
    