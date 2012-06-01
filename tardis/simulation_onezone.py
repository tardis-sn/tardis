import plasma
import initialize
import atomic
import line
import photon
import constants
import numpy as np
import montecarlo_multizone
import time
import model

def run_onezone(conn, fname):
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


def run_onezone_macro(conn, fname):
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

def run_multizone(conn, fname, max_atom=30, max_ion=30):
    initial_config = initialize.read_simple_config(fname)

    
    
    
    max_atom = 30
    max_ion = 30
    
    abundances = plasma.named2array_abundances(initial_config['abundances'], max_atom)
    atomic_model = atomic.KuruczMacroAtomModel.from_db(conn)
    line_list = line.read_line_list(conn)
    line2level = line_list['global_level_id_upper'] - 1
    atomic_model.macro_atom.read_nus(line_list['nu'])
    
    w7model = model.MultiZoneRadial.from_w7(initial_config['time_exp'])
    t_rad = initial_config['t_outer']
    t_inner = t_rad
    surface_inner = 4 * np.pi * w7model.r_inner[0]**2
    w = 1.
    energy_of_packet = 1. / initial_config['packets']
    volume = (4./3) * np.pi * (w7model.r_outer**3 - w7model.r_inner**3)
    w7model.set_line_list(line_list)
    w7model.set_atomic_model(atomic_model)
    w7model.read_abundances_uniform(abundances)
    w7model.initialize_plasmas(t_rad)
    
    i = 0
    print "starting", t_rad
    track_ws = []
    track_t_rads = []
    track_t_inner = []
    while True:
        start_time = time.time()
        track_t_rads.append(w7model.t_rads.copy())
        track_t_inner.append(t_inner)
        track_ws.append(w7model.ws.copy())
        i+=1
        if i > 5: break
        
        #return sn_plasma
        
        tau_sobolevs = w7model.calculate_tau_sobolevs()
        transition_probabilities = w7model.calculate_transition_probabilities(tau_sobolevs)

        nu_input = photon.random_blackbody_nu(t_inner, nu_range=(1e8*constants.c / 1, 1e8*constants.c/ 100000.), size=initial_config['packets'])
        mu_input = np.sqrt(np.random.random(initial_config['packets']))

        print "calculating transition probabilities"
        
        
        print "running montecarlo"
        nu, energy, nu_reabsorbed, energy_reabsorbed, j_estimators, nubar_estimators = \
                                    montecarlo_multizone.run_simple_oned(nu_input,
                                                                mu_input,
                                                                line_list['nu'],
                                                                tau_sobolevs,
                                                                w7model.r_inner,
                                                                w7model.r_outer,
                                                                w7model.v_inner,
                                                                w7model.electron_densities,
                                                                energy_of_packet,
                                                                transition_probabilities,
                                                                w7model.atomic_model.macro_atom.transition_type_total,
                                                                w7model.atomic_model.macro_atom.target_level_total,
                                                                w7model.atomic_model.macro_atom.target_line_total,
                                                                w7model.atomic_model.macro_atom.level_references,
                                                                line2level)
        
        new_t_rads = constants.trad_estimator_constant * nubar_estimators / j_estimators
        new_ws = constants.w_estimator_constant * j_estimators / (initial_config['time_simulation'] * new_t_rads**4 * volume)
        
        emitted_energy_fraction = np.sum(energy[nu != 0]) / 1.
        
        
        new_t_inner = (initial_config['luminosity_outer'] / (emitted_energy_fraction * constants.sigma_sb * surface_inner ))**.25
        
        for new_t_rad, new_w, old_trad, old_w in zip(new_t_rads, new_ws, w7model.t_rads, w7model.ws):
            print "new_trad / t_rad %.2f / %.2f, new_w/w %.5f %.5f" % (new_t_rad, old_trad, new_w, old_w)
        print '\n\n---'
        print "new_tinner / t_inner  %.2f / %.2f" %(new_t_inner, t_inner)
        w7model.ws = (new_ws + w7model.ws) * .5
        w7model.t_rads = (new_t_rads + w7model.t_rads) * .5
        t_inner = (new_t_inner + t_inner) * .5
        w7model.update_model(w7model.t_rads, w7model.ws)
        
        print "took %.2f s" % (time.time() - start_time)
        
    return nu_input, energy_of_packet, nu, energy, nu_reabsorbed, energy_reabsorbed, track_t_rads, track_ws, track_t_inner
    