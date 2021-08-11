from numba import njit
from tardis.montecarlo.montecarlo_numba import njit_dict, njit_dict_no_parallel
from tardis.montecarlo.montecarlo_numba.numba_interface import (
    LineInteractionType,
)

from tardis.montecarlo import (
    montecarlo_configuration as montecarlo_configuration,
)
from tardis.montecarlo.montecarlo_numba.frame_transformations import (
    get_doppler_factor,
    get_inverse_doppler_factor,
    angle_aberration_CMF_to_LF,
)
from tardis.montecarlo.montecarlo_numba.r_packet import (
    InteractionType,
)
from tardis.montecarlo.montecarlo_numba.utils import get_random_mu
from tardis.montecarlo.montecarlo_numba.macro_atom import (
        macro_atom, MacroAtomTransitionType
)



def scatter(r_packet, time_explosion):

    old_doppler_factor = get_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion
    )
    comov_nu = r_packet.nu * old_doppler_factor
    comov_energy = r_packet.energy * old_doppler_factor
    r_packet.mu = get_random_mu()
    inverse_new_doppler_factor = get_inverse_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion
    )

    r_packet.energy = comov_energy * inverse_new_doppler_factor
 
    return comov_nu, inverse_new_doppler_factor

# Maybe make the continuum selection a method of the continuum?
@njit(**njit_dict_no_parallel)
def continuum_event(r_packet, time_explosion, continuum, numba_plasma):

    comov_nu, inverse_new_doppler_factor = scatter(r_packet, time_explosion)

    old_doppler_factor = get_doppler_factor(r_packet.r, r_packet.mu, time_explosion)
    comov_nu = r_packet.nu * old_doppler_factor

    zrand = np.random.random()

    # Does this need to be re-calculated?
    continuum.calculate(comov_nu, r_packet.current_shell_id)

    # Need to determine if a collisional process or not.  If not:

    # This somehow determines the transition type and continuum id needed
    # but does not sample new frequencies
    destination_level_idx = continuum.determine_macro_activation_idx(
            comov_nu, r_packet.current_shell_id)

    transition_id, transition_type = macro_atom(
            destination_level_idx, 
            r_packet.current_shell_id, 
            numba_plasma
            )

    # Then use this to get a transition id from the macroatom 
    if transition_type == MacroAtomTransitionType.FF_EMISSION: 
        free_free_emission(r_packet, time_explosion, numba_plasma)
    elif transition_type = MacroAtomTransitionType.BF_EMISSION:
        bound_free_emission(r_packet, 
                time_explosion, 
                numba_plasma, 
                continuum, 
                transition_id)

# TODO: numbafy | Add cooling rates to numba plasma
def get_emission_probabilities(plasma, shell):
    
    C_fb = sim.plasma.cool_rate_fb_tot.iloc[0, shell]
    C_ff = sim.plasma.cool_rate_ff.iloc[0, shell]
    C_cl = sim.plasma.cool_rate_adiabatic.iloc[0, shell]
    C_cl += sim.plasma.cool_rate_coll_ion.iloc[0, shell]
    
    pi_fb = C_fb / (C_fb + C_ff + C_cl)
    pi_ff = C_ff / (C_fb + C_ff + C_cl)
    
    return pi_fb, pi_ff

# TODO: numbafy
def free_free_absorption(r_packet, numba_plasma):
    # Do the k_packet thing, but don't actually make one
    
    pi_fb, pi_ff = get_emission_probabilities(numba_plasma, shell)
    
    # a little hack
    z = np.random.random() - pi_fb
    if z < 0:
        free_bound_emission(r_packet, time_explosion, numba_plasma) 
    elif z < pi_ff:
        free_free_emission(r_packet, time_explosion, numba_plasma)
    else:
        macroatom()

@njit(**njit_dict_no_parallel)
def get_current_line_id(nu, numba_plasma):
    '''
    Get the next line id corresponding to a packet at frequency nu
    '''

    reverse_line_list = numba_plasma.line_list_nu[::-1]
    number_of_lines = len(numba_plasma.line_list_nu)
    line_id = number_of_lines - np.searchsorted(reverse_line_list, nu)
    return line_id


@njit(**njit_dict_no_parallel)
def free_free_emission(r_packet, time_explosion, numba_plasma):
    
    inverse_doppler_factor = get_inverse_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion
    )
    # Need to get the sampler into numba somehow
    # maybe I'll update the numba_plasma?
    comov_nu = numba_plasma.nu_ff_sampler(r_packet.current_shell_id)
    r_packet.nu = comov_nu * inverse_doppler_factor
    current_line_id = get_current_line_id(r_packet.nu, numba_plasma) 
    r_packet.next_line_id = current_line_id
    
    if montecarlo_configuration.full_relativity:
        r_packet.mu = angle_aberration_CMF_to_LF(
            r_packet, time_explosion, r_packet.mu
        )

@njit(**njit_dict_no_parallel)
def free_bound_emission(r_packet, time_explosion, numba_plasma, continuum, continuum_d):

    inverse_doppler_factor = get_inverse_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion
    )

    comov_nu = numba_plasma.nu_fb_sampler(
            r_packet.current_shell_id, 
            continuum_id
            )
    
    r_packet.nu = comov_nu * inverse_doppler_factor
    current_line_id = get_current_line_id(r_packet.nu, numba_plasma)
    r_packet.next_line_id = current_line_id

    if montecarlo_configuration.full_relativity:
        r_packet.mu = angle_aberration_CMF_to_LF(
            r_packet, time_explosion, r_packet.mu
        )



@njit(**njit_dict_no_parallel)
def thomson_scatter(r_packet, time_explosion):
    """
    Thomson scattering — no longer line scattering
    \n1) get the doppler factor at that position with the old angle
    \n2) convert the current energy and nu into the comoving frame with the old mu
    \n3) Scatter and draw new mu - update mu
    \n4) Transform the comoving energy and nu back using the new mu

    Parameters
    ----------
    r_packet : tardis.montecarlo.montecarlo_numba.r_packet.RPacket
    time_explosion : float
        time since explosion in seconds
    """
    old_doppler_factor = get_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion
    )
    comov_nu = r_packet.nu * old_doppler_factor
    comov_energy = r_packet.energy * old_doppler_factor
    r_packet.mu = get_random_mu()
    inverse_new_doppler_factor = get_inverse_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion
    )

    r_packet.nu = comov_nu * inverse_new_doppler_factor
    r_packet.energy = comov_energy * inverse_new_doppler_factor
    if montecarlo_configuration.full_relativity:
        r_packet.mu = angle_aberration_CMF_to_LF(
            r_packet, time_explosion, r_packet.mu
        )


@njit(**njit_dict_no_parallel)
def line_scatter(r_packet, time_explosion, line_interaction_type, numba_plasma):
    """
    Line scatter function that handles the scattering itself, including new angle drawn, and calculating nu out using macro atom

    Parameters
    ----------
    r_packet : tardis.montecarlo.montecarlo_numba.r_packet.RPacket
    time_explosion : float
    line_interaction_type : enum
    numba_plasma : tardis.montecarlo.montecarlo_numba.numba_interface.NumbaPlasma
    """

    old_doppler_factor = get_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion
    )
    r_packet.mu = get_random_mu()

    inverse_new_doppler_factor = get_inverse_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion
    )

    comov_energy = r_packet.energy * old_doppler_factor
    r_packet.energy = comov_energy * inverse_new_doppler_factor

    if line_interaction_type == LineInteractionType.SCATTER:
        line_emission(
            r_packet, r_packet.next_line_id, time_explosion, numba_plasma
        )
    else:  # includes both macro atom and downbranch - encoded in the transition probabilities
        activation_level_id = numba_plasma.line2macro_level_upper[
            r_packet.next_line_id
        ]
        emission_line_id, transition_type = macro_atom(
            activation_level_id,
            r_packet.current_shell_id,
            numba_plasma
        )
        line_emission(r_packet, emission_line_id, time_explosion, numba_plasma)


@njit(**njit_dict_no_parallel)
def line_emission(r_packet, emission_line_id, time_explosion, numba_plasma):
    """
    Sets the frequency of the RPacket properly given the emission channel

    Parameters
    ----------
    r_packet : tardis.montecarlo.montecarlo_numba.r_packet.RPacket
    emission_line_id : int
    time_explosion : float
    numba_plasma : tardis.montecarlo.montecarlo_numba.numba_interface.NumbaPlasma
    """
    r_packet.last_line_interaction_out_id = emission_line_id

    if emission_line_id != r_packet.next_line_id:
        pass
    inverse_doppler_factor = get_inverse_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion
    )
    r_packet.nu = (
        numba_plasma.line_list_nu[emission_line_id] * inverse_doppler_factor
    )
    r_packet.next_line_id = emission_line_id + 1
    nu_line = numba_plasma.line_list_nu[emission_line_id]

    if montecarlo_configuration.full_relativity:
        r_packet.mu = angle_aberration_CMF_to_LF(
            r_packet, time_explosion, r_packet.mu
        )
