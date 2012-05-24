# running the simulation
import photon
import model
import constants
import line
import numpy as np
import pdb
import operator

np.random.seed(25081980)

def runSimpleOneD(photons, mus, line_list, tau_sobolevs, model):
    atmosphere = model
    #line list: nu , oscillator strength, dens_line
    line_list_nu = constants.c / line_list
    nu = []
    energy = []
    scatter_events = 0
    e_scatter_events = 0
    resonance_events = 0
    reabsorb_counter =0
    photon_loop_counter = 0
    packet_no = len(photons)
    
    for i in xrange(packet_no):
        if i%100 == 0: print '@ packet %d' % i
        #print '@ packet %d' % i
        #copy by value!!
        current_r = atmosphere.r_inner
        reabsorbed = None
        current_nu = photons[i]
        current_mu = mus[i]
        
        
        
        current_energy = 1.
        
        
        current_nu_cmf = current_nu * (1. - (current_mu * model.v_inner * constants.inverse_c))
        cur_line_id = line_list_nu.size - line_list_nu[::-1].searchsorted(current_nu_cmf)
        if cur_line_id == line_list_nu.size: last_line=True
        else: last_line = False

        assert cur_line_id >= 0
        loop_counter = 0
        tau_event = -np.log(np.random.random())
        
        close_line = False
        interaction_history = []
        current_nu_history = []
        current_nu_cmf_history = []
        current_line_history = []
        current_mu_history = []
        current_r_history = []
        d_line_history = []
        #pdb.set_trace()
        while True:
            loop_counter +=1
            
            if not last_line:
                nu_line  = line_list_nu[cur_line_id]
            
            if close_line == True:
                d_line = 0.0
                interaction_type = 3
                close_line = False
            else:
                d_inner = atmosphere.compute_distance2inner(current_r, current_mu)
                d_outer = atmosphere.compute_distance2outer(current_r, current_mu)
                if last_line:
                    d_line = 1e99
                else:
                    d_line = atmosphere.compute_distance2line(current_r, current_mu, current_nu, nu_line)
                d_line_history.append(d_line)
                d_electron = atmosphere.compute_distance2electron(current_r, current_mu, tau_event)
                interaction_type = min(enumerate([d_inner, d_outer, d_electron, d_line]), key=operator.itemgetter(1))[0]
                

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
                e_scatter_events += 1
                current_r = line.get_r_sobolev(current_r, current_mu, d_electron)
                nu_electron = current_nu * (1 - (current_mu * current_r * atmosphere.inverse_t_exp * constants.inverse_c))
                energy_electron = current_energy * (1 - (current_mu * current_r * atmosphere.inverse_t_exp * constants.inverse_c))
                current_mu = 2*np.random.random() - 1
                current_nu = nu_electron / (1 - (current_mu * current_r * atmosphere.inverse_t_exp * constants.inverse_c))
                current_energy = energy_electron / (1 - (current_mu * current_r * atmosphere.inverse_t_exp * constants.inverse_c))
                tau_event = -np.log(np.random.random())
            
            elif interaction_type == 3:    
                #It has a chance to hit the line
                # Compute the t_sob of the line
                
                tau_sobolev = tau_sobolevs[cur_line_id]
                
                resonance_events += 1
                prev_r = current_r
                current_r = line.get_r_sobolev(current_r, current_mu, d_line)
                
                cur_line_id += 1
                
                #check for last line
                if cur_line_id == line_list_nu.size:
                        last_line = True
                #check for same line        
                if not last_line:
                    if abs(line_list_nu[cur_line_id] - nu_line) < 0.1:
                        close_line = True
                    
                #Check for line interaction
                if tau_event < tau_sobolev:
                    #print 'scatter'
                    #line event happens - move and scatter packet
                    #choose new mu
                    current_mu = 2*np.random.random() - 1
                    doppler_factor = 1/(1 - (current_mu * current_r * atmosphere.inverse_t_exp * constants.inverse_c))
                    current_nu = nu_line * doppler_factor
                    current_energy *= doppler_factor
                    #current_energy = 0
                    tau_event = -np.log(np.random.random())
                    scatter_events += 1
                    scattered = True
                    
                else:
                    scattered = False
                    prev_mu = current_mu
                    if d_line != 0.0:
                        current_mu = (d_line**2 + current_r**2 - prev_r**2) / (2*d_line*current_r)
                    tau_event -= tau_sobolev
                
            
        assert (reabsorbed == True or reabsorbed == False)
        
        if reabsorbed == True:
            reabsorb_counter += 1
            continue
        elif reabsorbed == False:
            nu.append(current_nu)
            energy.append(current_energy)
    print 'resonance, scatter, photon_loops', resonance_events, scatter_events, photon_loop_counter
    print '%d packets were absorbed' % reabsorb_counter
    return nu, energy