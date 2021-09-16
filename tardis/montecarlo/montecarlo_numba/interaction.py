from numba import njit
import numpy as np
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
    InteractionType, PacketStatus
)
from tardis.montecarlo.montecarlo_numba.utils import get_random_mu
from tardis.montecarlo.montecarlo_numba.macro_atom import (
        macro_atom, MacroAtomTransitionType
)



@njit(**njit_dict_no_parallel)
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


    # Does this need to be re-calculated?
    # continuum.calculate(comov_nu, r_packet.current_shell_id)

    # Need to determine if a collisional process or not.  If not:

    # This somehow determines the transition type and continuum id needed
    # but does not sample new frequencies
    # This reruns continuum.calculate

    #doppler_factor = get_doppler_factor(r_packet.r, r_packet.mu, time_explosion)
    #comov_nu = r_packet.nu * doppler_factor
    old_doppler_factor = get_doppler_factor(
            r_packet.r, 
            r_packet.mu, 
            time_explosion
            )

    #comov_nu = r_packet.nu * old_doppler_factor
    r_packet.mu = get_random_mu()
    inverse_doppler_factor = get_inverse_doppler_factor(
            r_packet.r, 
            r_packet.mu, 
            time_explosion
            )
    comov_energy = r_packet.energy * old_doppler_factor
    r_packet.energy = comov_energy * inverse_doppler_factor

    destination_level_idx = continuum.determine_macro_activation_idx(
            r_packet.nu, r_packet.current_shell_id)

    macro_atom_event(destination_level_idx, 
            r_packet, time_explosion,
            numba_plasma, continuum)

@njit(**njit_dict_no_parallel)
def macro_atom_event(destination_level_idx, 
    r_packet, time_explosion, numba_plasma, continuum):

    try:
        transition_id, transition_type = macro_atom(
                destination_level_idx, 
                r_packet.current_shell_id, 
                numba_plasma
                )
    except:
        print('destination_level', destination_level_idx)
        print('r_packet.nu', r_packet.nu)
        raise Exception("Macroatom ran out of block")

    # Then use this to get a transition id from the macroatom 
    if transition_type == MacroAtomTransitionType.FF_EMISSION: 
        print("FF EMISSION")
        free_free_emission(
                r_packet, 
                time_explosion, 
                numba_plasma, 
                continuum
                )
    
    elif transition_type == MacroAtomTransitionType.BF_EMISSION:
        print("BF EMISSION")
        bound_free_emission(
                r_packet, 
                time_explosion, 
                numba_plasma, 
                continuum, 
                transition_id
                )
    elif transition_type == MacroAtomTransitionType.BF_COOLING:
        print("BF COOLING")
        #print(' before r_packet.nu=', r_packet.nu, r_packet.next_line_id, numba_plasma.line_list_nu[r_packet.next_line_id])
        bf_cooling(r_packet, time_explosion, numba_plasma, continuum)
        #free_free_emission(
        #        r_packet, 
        #        time_explosion, 
        #        numba_plasma, 
        #        continuum
        #        )
        #print(' now r_packet.nu=', r_packet.nu)
    
    elif transition_type == MacroAtomTransitionType.ADIABATIC_COOLING:
        print("ADIABATIC")
        adiabatic_cooling(
                r_packet, 
                time_explosion, 
                numba_plasma, 
                continuum
                )

    elif transition_type == MacroAtomTransitionType.BB_EMISSION:
        print("LINE EMISSION")
        line_emission(
                r_packet, 
                transition_id,
                time_explosion,
                numba_plasma
                )
    else:
        print("OTHER CONTINUUM")
        pass


@njit(**njit_dict_no_parallel)
def bf_cooling(r_packet, time_explosion, numba_plasma, continuum):

    fb_cooling_prob = numba_plasma.p_fb_deactivation[:, r_packet.current_shell_id]
    #print('fb_cooling_prob length=', len(fb_cooling_prob))
    p = fb_cooling_prob[0]
    i = 0
    zrand = np.random.random()
    #print('zrand=',zrand)
    #print('sum=', fb_cooling_prob.sum())
    while p <= zrand:
        i += 1
        p += fb_cooling_prob[i]
    #print('Got continuum_idx=',i)
    continuum_idx = i
    bound_free_emission(
            r_packet,
            time_explosion,
            numba_plasma,
            continuum,
            continuum_idx
            )

@njit(**njit_dict_no_parallel)
def adiabatic_cooling(r_packet, time_explosion, numba_plasma, continuum):

    r_packet.status = PacketStatus.REABSORBED
    # This should update the heating estimators by addiing the comoving energy
    # Then go into a k-packet creation
    pass

@njit(**njit_dict_no_parallel)
def get_current_line_id(nu, line_list):
    '''
    Get the next line id corresponding to a packet at frequency nu
    '''

    reverse_line_list = line_list[::-1]
    number_of_lines = len(line_list)
    line_id = number_of_lines - np.searchsorted(reverse_line_list, nu)
    return line_id


@njit(**njit_dict_no_parallel)
def free_free_emission(r_packet, time_explosion, numba_plasma, continuum):
    #print("Free Free emission")
    inverse_doppler_factor = get_inverse_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion
    )
    comov_nu = continuum.sample_nu_free_free(r_packet.current_shell_id)
    r_packet.nu = comov_nu * inverse_doppler_factor
    current_line_id = get_current_line_id(
            comov_nu,
            numba_plasma.line_list_nu
            ) 
    r_packet.next_line_id = current_line_id
    
    if montecarlo_configuration.full_relativity:
        r_packet.mu = angle_aberration_CMF_to_LF(
            r_packet, time_explosion, r_packet.mu
        )
    #print("End Free Free emission")

@njit(**njit_dict_no_parallel)
def bound_free_emission(r_packet, time_explosion, numba_plasma, continuum, continuum_id):

    #print("sampling bf")
    #print('r_packet id = ', r_packet.index)
    #if not len(continuum.current_continua):
    #    #print("Error: No Current Continuum!")
    #    return
    #print("BF emission")
    #print(" r_packet.nu", r_packet.nu)
    inverse_doppler_factor = get_inverse_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion
    )

    comov_nu = continuum.sample_nu_free_bound(
            r_packet.current_shell_id, 
            continuum_id)
    #print('  Sampled nu_fb, got:', comov_nu)
    r_packet.nu = comov_nu * inverse_doppler_factor
    current_line_id = get_current_line_id(
            comov_nu, 
            numba_plasma.line_list_nu
            )
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
    #print("Thomson scatter")
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

    #print("End Thomson scatter")

@njit(**njit_dict_no_parallel)
def line_scatter(r_packet, time_explosion, 
        line_interaction_type, numba_plasma, continuum):
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
        macro_atom_event(activation_level_id, r_packet, 
                time_explosion, numba_plasma, continuum)
        #emission_line_id, transition_type = macro_atom(
        #    activation_level_id,
        #    r_packet.current_shell_id,
        #    numba_plasma
        #)
        #line_emission(r_packet, emission_line_id, time_explosion, numba_plasma)


@njit(**njit_dict_no_parallel)
def line_emission(r_packet, emission_line_id, time_explosion, 
        numba_plasma):
    """
    Sets the frequency of the RPacket properly given the emission channel

    Parameters
    ----------
    r_packet : tardis.montecarlo.montecarlo_numba.r_packet.RPacket
    emission_line_id : int
    time_explosion : float
    numba_plasma : tardis.montecarlo.montecarlo_numba.numba_interface.NumbaPlasma
    """
    #print("BB emission")
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
    #print("end BB emission")
