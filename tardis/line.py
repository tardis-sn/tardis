#calculation for line interaction
import constants
import numpy as np
def get_tau_line(nu, f_line, dens_line, t_exp):
    tau_line = constants.sobolev_coeff * f_line * (constants.c / nu) * dens_line * t_exp
    return tau_line

def get_r_sobolev(r, mu, d_line):
    return np.sqrt(r**2 + d_line**2 + 2 * r * d_line * mu)