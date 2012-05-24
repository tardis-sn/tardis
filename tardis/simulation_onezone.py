import plasma
import initialize
import line
def run_oned(conn, fname):
    initial_config = initialize.read_simple_config(fname)
    return initial_config
    
    t_rad = initial_config['t_inner']
    w = 1.
    energy_of_packet = 1 / initial_config['packets']
    
    sn_plasma = plasma.NebularPlasma.from_db(initial_config['abundances'], inital_config['density'], conn)
    
    line_list = line.read_line_list(conn)
    
    while True:
        sn_plasma.update_radiationfield(t_rad = t_rad)
        tau_sobolev = sn_plasma.calculate_tau_sobolev(line_list, initial_config['time_exp'])
        
        nus = photon.random_blackbody_nu(T, nu_range=(1e8*constants.c / 2000, 1e8*constants.c/ 9000), size=no_of_packets)
        mus = np.sqrt(np.random.random(no_of_packets))
        
        #WRITE!!
        tau_sobolev = myplasma.calculate_tau_sobolev()
        
        
        #add energy to montecarlo function
        nu, energy, j_estimator, nubar_estimator = c_montecarlo.run_simple_oned(nus, mus, line_list_nu, tau_sobolevs,
                    initial_config['r_inner'], initial_config['r_outer'], initial_config['v_inner'], myplasma.electron_density, energy_of_packet)
        
        new_trad = constants.trad_estimator_constant * nubar_estimator / j_estimator
        new_w = constants.w_estimator_constant * j_estimator / new_trad**4 / (inital_config['time_simulation'] * (4/3) * np.pi * (r_outer**3 - r_inner**3))
        print "trad / new_trad = %.2f/%.2f w / new_w =  %.2f/%.2f" % (trad, new_trad, w, new_w)
        w = (new_w + w)*.5
        trad = (new_trad + trad)*.5
        
        