# running the simulation
import photon
import model
import constants
import numpy as np

def runSimpleOneD(packets=1000000, phi=1e-13):
    nu_line = constants.c / 1216.10 / 1e-8;
    #nu_line = constants.c / 1500.10 / 1e-8;
    atmosphere = model.OneZoneBaseModel(6.96e10, 10*6.96e10, 0.01*constants.c)
    photon_source = photon.SimplePhotonSource.from_wavelength(1000, 2000)
    nu = []
    energy = []
    for i in xrange(packets):
        if i%10000 == 0: print '@ packet %d' % i
        #print '@ packet %d' % i
        #copy by value!!
        current_r = atmosphere.r_inner
        reabsorbed = None
        current_nu, current_mu = photon_source()
        current_energy = 1.
        while True:
            d_inner = atmosphere.compute_distance2inner(current_r, current_mu)
            d_outer = atmosphere.compute_distance2outer(current_r, current_mu)
            d_line = atmosphere.compute_distance2line(current_r, current_mu, current_nu, nu_line)
            
            #packet escaping
            if (d_outer < d_line) and (d_outer < d_inner):
                reabsorbed = False
                break
            
            #packet reabosrbing into core
            elif (d_inner < d_line) and (d_inner < d_outer):
                reabsorbed = True
                break
            
            #line scattering event
            elif (d_line < d_inner) and (d_line < d_outer):
                #It has a chance to hit the line
                # Compute the t_sob of the line
                r_sobolev = np.sqrt( current_r**2 + d_line**2 + 2 * current_r * d_line * current_mu)
                #nh_sob = phi /4 /PI/r_sob/r_sob/r_sob/v1*r1/MH*MSUN/YEAR;
                nh_sob = (phi / (4 * 3.1 * r_sobolev**3 * atmosphere.v_outer)) * (atmosphere.r_outer * 1.9e33) / (1.67e-24 * 3.15e7)
                #nh_sob = (1e-13 / (4 * 3.1 * r_sobolev**3 * atmosphere.v_outer)) * (atmosphere.r_outer * 1.9e33) / (1.67e-24 * 3.15e7)
                #nh_sob = 0
                tau_sobolev = 0.02655103 * 0.416 * nh_sob * constants.c * atmosphere.r_outer / current_nu / atmosphere.v_outer
                tau_random = -np.log(np.random.random())
                #Check for line interaction
                if tau_random > tau_sobolev:
                # No line event check if it gets reabsorbed or leaves
                    if d_inner > d_outer:
                        reabsorbed = False
                        break
                    else:
                        reabsorbed = True
                        break
                else:
                    #line event happens - move and scatter packet
                    current_r = r_sobolev
                    #choose new mu
                    current_mu = 2*np.random.random() - 1
                    current_nu = nu_line / (1 - (current_mu * current_r * atmosphere.homol_coeff / constants.c))
                    current_energy /= (1 - (current_mu * current_r * atmosphere.homol_coeff / constants.c))
        assert (reabsorbed == True or reabsorbed == False)
        if reabsorbed == True:
            continue
        elif reabsorbed == False:
            nu.append(current_nu)
            energy.append(current_energy)
    return nu, energy