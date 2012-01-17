# running the simulation
import photon
import model
import constants
import line
import numpy as np
#import pdb
import operator

np.random.seed(25081980)
def runSimpleOneD(packets=1000000, phi=1e-13):
    atmosphere = model.OneZoneBaseModel(6.96e10, 10*6.96e10, 0.01*constants.c)
    photon_source = photon.SimplePhotonSource.from_wavelength(1000, 2000)
    #line list: nu , oscillator strength, dens_line
    line_list = [(constants.c / 1100. / 1e-8, 1, 1e10), 
        (constants.c / 1200.10 / 1e-8, 1, 1e10)]
    line_list_nu = np.array(zip(*line_list[:-1])[0][::-1])
    nu = []
    energy = []
    scatter_events = 0
    resonance_events = 0
    photon_loop_counter = 0
    wl = np.arange(1000, 2000,1)
    
    for i in xrange(packets):
        
        if i%1000 == 0: print '@ packet %d' % i
        #print '@ packet %d' % i
        #copy by value!!
        current_r = atmosphere.r_inner
        reabsorbed = None
        current_nu, current_mu = photon_source()
        #current_nu = constants.c / wl[i] / 1e-8
        #current_mu = 1
        current_nu_cmf = current_nu * (1. - (current_mu * atmosphere.v_inner * constants.inverse_c))
        current_energy = 1.
        cur_line_id = line_list_nu.size - line_list_nu.searchsorted(current_nu_cmf)
        #print "wl, cur_line_id", 1e8 * constants.c / current_nu, cur_line_id,
        assert cur_line_id >=0
        loop_counter = 0
        tau_event = -np.log(np.random.random())
        last_line = False
        #pdb.set_trace()
        while True:
            loop_counter +=1
            nu_line, f_line, dens_line = line_list[cur_line_id]
            d_inner = atmosphere.compute_distance2inner(current_r, current_mu)
            d_outer = atmosphere.compute_distance2outer(current_r, current_mu)
            if last_line:
                d_line = 1e99
            else:
                d_line = atmosphere.compute_distance2line(current_r, current_mu, current_nu, nu_line)
            
            d_electron = atmosphere.compute_distance2electron(current_r, current_mu, tau_event)
            
            interaction_type = min(enumerate([d_inner, d_outer, d_electron, d_line]), key=operator.itemgetter(1))[0]
            #packet escaping
            assert abs(d_line > 1.)
            #if (d_outer < d_inner) and (d_outer < d_electron) and (d_outer < d_line):
            if interaction_type == 0:
                #reabsorbed
                reabsorbed = True
                break
            
            #packet reabosrbing into core
            elif interaction_type == 1:
            #elif (d_inner < d_outer) and (d_inner < d_electron) and (d_inner < d_line):
                #escaped
                reabsorbed = False
                break
            elif interaction_type == 2:
            #elif (d_electron < d_outer) and (d_electron < d_inner) and (d_electron < d_line):
            #electron scattering
                
                current_r = line.get_r_sobolev(current_r, current_mu, d_electron)
                nu_electron = current_nu * (1 - (current_mu * current_r * atmosphere.inverse_t_exp * constants.inverse_c))
                energy_electron = current_energy * (1 - (current_mu * current_r * atmosphere.inverse_t_exp * constants.inverse_c))
                current_mu = 2*np.random.random() - 1
                current_nu = nu_electron / (1 - (current_mu * current_r * atmosphere.inverse_t_exp * constants.inverse_c))
                current_energy = energy_electron / (1 - (current_mu * current_r * atmosphere.inverse_t_exp * constants.inverse_c))
                tau_event = -np.log(np.random.random())
            
            elif interaction_type == 3:    
            #elif (d_line < d_outer) and (d_line < d_electron) and (d_line < d_inner):
                #print "resonance", cur_line_id, 
                #It has a chance to hit the line
                # Compute the t_sob of the line
                tau_sobolev = line.get_tau_line(nu_line, f_line, dens_line, atmosphere.t_exp)
                #print tau_sobolev
                current_r = line.get_r_sobolev(current_r, current_mu, d_line)
                if cur_line_id == (len(line_list) -1):
                    last_line = True
                else:
                    cur_line_id +=1
                #Check for line interaction
                if tau_event < tau_sobolev:
                    #print 'scatter'
                    #line event happens - move and scatter packet
                    #choose new mu
                    current_mu = 2*np.random.random() - 1
                    current_nu = nu_line / (1 - (current_mu * current_r * atmosphere.inverse_t_exp * constants.inverse_c))
                    current_energy /= (1 - (current_mu * current_r * atmosphere.inverse_t_exp * constants.inverse_c))
                    tau_event = -np.log(np.random.random())
                    scatter_events += 1
                else:
                    tau_event -= tau_sobolev
            photon_loop_counter = max(photon_loop_counter, loop_counter)
        assert (reabsorbed == True or reabsorbed == False)
        if reabsorbed == True:
            continue
        elif reabsorbed == False:
            nu.append(current_nu)
            energy.append(current_energy)
    print 'resonance, scatter, photon_loops', resonance_events, scatter_events, photon_loop_counter
    return nu, energy