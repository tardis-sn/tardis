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
from tardis.montecarlo.montecarlo_numba.macro_atom import macro_atom

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

def continuum_event(r_packet, time_explosion, continuum):

    comov_nu, inverse_new_doppler_factor = scatter(r_packet, time_explosion)

    old_doppler_factor = get_doppler_factor(r_packet.r, r_packet.mu, time_explosion)
    comov_nu = r_packet.nu * old_doppler_factor

    zrand = np.random.random()

    # Does this need to be re-calculated?
    continuum.calculate(comov_nu, r_packet.current_shell_id)
    chi_continuum = continuum.chi_bf + continuum.chi_ff

    # Since trace_packet differentiates between thomson scattering
    # and other continuum processes, we need to renormalize
    # our odds of selecting each continuum process given that
    # we are not thomson scattering
    # P(bf|~e_scat) = P(~e_scat|bf) P(bf) / P(~e_scat)
    # P(~e_scat) = 1 - P(e_scat) = 1 - chi_e / chi_nu
    # P(bf) = chi_bf / chi_nu
    # P(bf|~e_scat) = chi_bf / chi_nu / (1 - chi_e / chi_nu)
    # P(bf|~e_scat) = chi_bf / (chi_nu - chi_e)
    # Since chi_nu is the sum of ff, bf, and e_scat opacities
    # we can just set chi_continuum = chi_nu - chi_e = chi_bf + chi_ff
    # alternatively, this could again just be selected in
    # trace packet, it seems pretty straightforward to do

    if zrand < chi_bf / chi_continuum:
        bound_free_absorption(r_packet, time_explosion, plasma)
    else:
        free_free_absorption(r_packet, time_explosion)


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
