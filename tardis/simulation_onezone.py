import plasma
import initialize
import atomic
import line
import photon
import constants
import numpy as np
import c_montecarlo
import time

def run_oned(conn, fname):
    initial_config = initialize.read_simple_config(fname)
    #return initial_config
    
    t_rad = initial_config['t_outer']
    t_inner = t_rad
    surface_inner = 4 * np.pi * initial_config['r_inner']**2
    w = 1.
    energy_of_packet = 1. / initial_config['packets']
    atomic_model = atomic.KuruczMacroAtomModel.from_db(conn)
    
    sn_plasma = plasma.NebularPlasma.from_model(initial_config['abundances'], initial_config['density'], atomic_model)
    
    line_list = line.read_line_list(conn, atoms=initial_config['abundances'].keys())
    
    print "%d lines used in simulation" % len(line_list)
    line_list_nu = constants.c / (line_list['wl'] * 1e-8)
    print "initial radiation temperature", t_rad
    volume = (4./3) * np.pi * (initial_config['r_outer']**3 - initial_config['r_inner']**3)
    i = 0
    while True:
        i+=1
        if i > 1: break
        sn_plasma.update_radiationfield(t_rad=t_rad, w=w)
        #return sn_plasma
        tau_sobolevs = sn_plasma.calculate_tau_sobolev(line_list, initial_config['time_exp'])
        
        nu_input = photon.random_blackbody_nu(t_inner, nu_range=(1e8*constants.c / 1, 1e8*constants.c/ 100000.), size=initial_config['packets'])
        mu_input = np.sqrt(np.random.random(initial_config['packets']))
        
        #WRITE!!
       
        
        
        #add energy to montecarlo function
        nu, energy, nu_reabsorbed, energy_reabsorbed, j_estimator, nubar_estimator = c_montecarlo.run_simple_oned(nu_input, mu_input, line_list_nu, tau_sobolevs,
                    initial_config['r_inner'], initial_config['r_outer'], initial_config['v_inner'], sn_plasma.electron_density, energy_of_packet)
        
        new_t_rad = constants.trad_estimator_constant * nubar_estimator / j_estimator
        new_w = constants.w_estimator_constant * j_estimator / (initial_config['time_simulation'] * new_t_rad**4 * volume)
        
        emitted_energy_fraction = np.sum(energy[nu != 0]) / 1.
        
        
        new_t_inner = (initial_config['luminosity_outer'] / (emitted_energy_fraction * constants.sigma_sb * surface_inner ))**.25
        
        print "trad / new_trad = %.2f / %.2f\n t_inner / new_t_inner %.2f / %.2f\n w / new_w =  %.5f/%.5f" %\
                                 (t_rad, new_t_rad, t_inner, new_t_inner, w, new_w)
        w = (new_w + w) * .5
        t_rad = (new_t_rad + t_rad) * .5
        t_inner = (new_t_inner + t_inner) * .5
        
    return nu_input, energy_of_packet, nu, energy, nu_reabsorbed, energy_reabsorbed


def run_oned_macro(conn, fname):
    initial_config = initialize.read_simple_config(fname)
    #return initial_config
    
    t_rad = initial_config['t_outer']
    t_inner = t_rad
    surface_inner = 4 * np.pi * initial_config['r_inner']**2
    w = 1.
    energy_of_packet = 1. / initial_config['packets']
    atomic_model = atomic.KuruczMacroAtomModel.from_db(conn)
    
    
    sn_plasma = plasma.NebularPlasma.from_model(initial_config['abundances'], initial_config['density'], atomic_model)
    
    line_list = line.read_line_list(conn)
    
    line2level = line_list['global_level_id_upper'] - 1
    
    
    print "%d lines used in simulation" % len(line_list)
    line_list_nu = constants.c / (line_list['wl'] * 1e-8)
    
    print "Initializing Macro atom probabilities"
    atomic_model.macro_atom.read_nus(line_list_nu)
    print "initial radiation temperature", t_rad
    volume = (4./3) * np.pi * (initial_config['r_outer']**3 - initial_config['r_inner']**3)
    i = 0
    while True:
        start_time = time.time()
        i+=1
        if i > 1: break
        sn_plasma.update_radiationfield(t_rad=t_rad, w=w)
        #return sn_plasma
        tau_sobolevs = sn_plasma.calculate_tau_sobolev(line_list, initial_config['time_exp'])
        
        
        print "calculating transition probabilities"
        p_transition = atomic_model.macro_atom.calculate_transition_probabilities(tau_sobolevs, t_rad, w)
        
        nu_input = photon.random_blackbody_nu(t_inner, nu_range=(1e8*constants.c / 1, 1e8*constants.c/ 100000.), size=initial_config['packets'])
        mu_input = np.sqrt(np.random.random(initial_config['packets']))
        
        #WRITE!!
       
        
        
        #add energy to montecarlo function
        print "running montecarlo"
        nu, energy, nu_reabsorbed, energy_reabsorbed, j_estimator, nubar_estimator = \
                                    c_montecarlo.run_simple_oned(nu_input,
                                                                mu_input,
                                                                line_list_nu,
                                                                tau_sobolevs,
                                                                initial_config['r_inner'],
                                                                initial_config['r_outer'],
                                                                initial_config['v_inner'],
                                                                sn_plasma.electron_density,
                                                                energy_of_packet,
                                                                p_transition,
                                                                atomic_model.macro_atom.transition_type_total,
                                                                atomic_model.macro_atom.target_level_total,
                                                                atomic_model.macro_atom.target_line_total,
                                                                atomic_model.macro_atom.level_references,
                                                                line2level)
        
        new_t_rad = constants.trad_estimator_constant * nubar_estimator / j_estimator
        new_w = constants.w_estimator_constant * j_estimator / (initial_config['time_simulation'] * new_t_rad**4 * volume)
        
        emitted_energy_fraction = np.sum(energy[nu != 0]) / 1.
        
        
        new_t_inner = (initial_config['luminosity_outer'] / (emitted_energy_fraction * constants.sigma_sb * surface_inner ))**.25
        
        print "trad / new_trad = %.2f / %.2f\n t_inner / new_t_inner %.2f / %.2f\n w / new_w =  %.5f/%.5f" %\
                                 (t_rad, new_t_rad, t_inner, new_t_inner, w, new_w)
        w = (new_w + w) * .5
        t_rad = (new_t_rad + t_rad) * .5
        t_inner = (new_t_inner + t_inner) * .5
        print "took %.2f s" % (time.time() - start_time)
    return nu_input, energy_of_packet, nu, energy, nu_reabsorbed, energy_reabsorbed