from numba import njit
import numpy as np
from tardis.montecarlo.montecarlo_numba import njit_dict, njit_dict_no_parallel
from tardis.montecarlo.montecarlo_numba.numba_interface import (
    LineInteractionType,
)

from tardis.montecarlo import (
    montecarlo_configuration as montecarlo_configuration,
)
from tardis.transport.frame_transformations import (
    get_doppler_factor,
    get_inverse_doppler_factor,
    angle_aberration_CMF_to_LF,
)
from tardis.montecarlo.montecarlo_numba.r_packet import (
    InteractionType,
    PacketStatus,
)
from tardis.montecarlo.montecarlo_numba.utils import get_random_mu
from tardis.montecarlo.montecarlo_numba.macro_atom import (
    macro_atom,
    MacroAtomTransitionType,
)
from tardis import constants as const

K_B = const.k_B.cgs.value
H = const.h.cgs.value


@njit(**njit_dict_no_parallel)
def determine_bf_macro_activation_idx(
    numba_plasma, nu, chi_bf_contributions, active_continua
):
    """
    Determine the macro atom activation level after bound-free absorption.

    Parameters
    ----------
    nu : float
        Comoving frequency of the r-packet.
    chi_bf_contributions : numpy.ndarray, dtype float
        Cumulative distribution of bound-free opacities at frequency
        `nu`.
    active_continua : numpy.ndarray, dtype int
        Continuum ids for which absorption is possible for frequency `nu`.

    Returns
    -------
    float
        Macro atom activation idx.
    """
    # Perform a MC experiment to determine the continuum for absorption
    index = np.searchsorted(chi_bf_contributions, np.random.random())
    continuum_id = active_continua[index]

    # Perform a MC experiment to determine whether thermal or
    # ionization energy is created
    nu_threshold = numba_plasma.photo_ion_nu_threshold_mins[continuum_id]
    fraction_ionization = nu_threshold / nu
    if (
        np.random.random() < fraction_ionization
    ):  # Create ionization energy (i-packet)
        destination_level_idx = numba_plasma.photo_ion_activation_idx[
            continuum_id
        ]
    else:  # Create thermal energy (k-packet)
        destination_level_idx = numba_plasma.k_packet_idx
    return destination_level_idx


@njit(**njit_dict_no_parallel)
def determine_continuum_macro_activation_idx(
    numba_plasma, nu, chi_bf, chi_ff, chi_bf_contributions, active_continua
):
    """
    Determine the macro atom activation level after a continuum absorption.

    Parameters
    ----------
    nu : float
        Comoving frequency of the r-packet.
    chi_bf : numpy.ndarray, dtype float
        Bound-free opacity.
    chi_bf : numpy.ndarray, dtype float
        Free-free opacity.
    chi_bf_contributions : numpy.ndarray, dtype float
        Cumulative distribution of bound-free opacities at frequency
        `nu`.
    active_continua : numpy.ndarray, dtype int
        Continuum ids for which absorption is possible for frequency `nu`.

    Returns
    -------
    float
        Macro atom activation idx.
    """
    fraction_bf = chi_bf / (chi_bf + chi_ff)
    # TODO: In principle, we can also decide here whether a Thomson
    # scattering event happens and need one less RNG call.
    if np.random.random() < fraction_bf:  # Bound-free absorption
        destination_level_idx = determine_bf_macro_activation_idx(
            numba_plasma, nu, chi_bf_contributions, active_continua
        )
    else:  # Free-free absorption (i.e. k-packet creation)
        destination_level_idx = numba_plasma.k_packet_idx
    return destination_level_idx


@njit(**njit_dict_no_parallel)
def sample_nu_free_free(numba_plasma, shell):
    """
    Attributes
    ----------
    nu_ff_sampler : float
        Frequency of the free-free emission process

    """
    T = numba_plasma.t_electrons[shell]
    zrand = np.random.random()
    return -K_B * T / H * np.log(zrand)


@njit(**njit_dict_no_parallel)
def sample_nu_free_bound(numba_plasma, shell, continuum_id):
    """
    Attributes
    ----------
    nu_fb_sampler : float
        Frequency of the free-bounds emission process
    """

    start = numba_plasma.photo_ion_block_references[continuum_id]
    end = numba_plasma.photo_ion_block_references[continuum_id + 1]
    phot_nus_block = numba_plasma.phot_nus[start:end]
    em = numba_plasma.emissivities[start:end, shell]

    zrand = np.random.random()
    idx = np.searchsorted(em, zrand, side="right")

    return phot_nus_block[idx] - (em[idx] - zrand) / (em[idx] - em[idx - 1]) * (
        phot_nus_block[idx] - phot_nus_block[idx - 1]
    )


@njit(**njit_dict_no_parallel)
def continuum_event(
    r_packet,
    time_explosion,
    numba_plasma,
    chi_bf_tot,
    chi_ff,
    chi_bf_contributions,
    current_continua,
):
    """
    continuum event handler - activate the macroatom and run the handler

    Parameters
    ----------
    r_packet : tardis.montecarlo.montecarlo_numba.r_packet.RPacket
    time_explosion : float
    numba_plasma : tardis.montecarlo.montecarlo_numba.numba_interface.NumbaPlasma
    continuum : tardis.montecarlo.montecarlo_numba.numba_interface.Continuum
    """
    old_doppler_factor = get_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion
    )

    r_packet.mu = get_random_mu()
    inverse_doppler_factor = get_inverse_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion
    )
    comov_energy = r_packet.energy * old_doppler_factor
    comov_nu = (
        r_packet.nu * old_doppler_factor
    )  # make sure frequency should be updated
    r_packet.energy = comov_energy * inverse_doppler_factor
    r_packet.nu = comov_nu * inverse_doppler_factor

    destination_level_idx = determine_continuum_macro_activation_idx(
        numba_plasma,
        comov_nu,
        chi_bf_tot,
        chi_ff,
        chi_bf_contributions,
        current_continua,
    )

    macro_atom_event(
        destination_level_idx, r_packet, time_explosion, numba_plasma
    )


@njit(**njit_dict_no_parallel)
def macro_atom_event(
    destination_level_idx, r_packet, time_explosion, numba_plasma
):
    """
    Macroatom event handler - run the macroatom and handle the result

    Parameters
    ----------
    destination_level_idx : int
    r_packet : tardis.montecarlo.montecarlo_numba.r_packet.RPacket
    time_explosion : float
    numba_plasma : tardis.montecarlo.montecarlo_numba.numba_interface.NumbaPlasma
    """

    transition_id, transition_type = macro_atom(
        destination_level_idx, r_packet.current_shell_id, numba_plasma
    )

    if (
        montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED
        and transition_type == MacroAtomTransitionType.FF_EMISSION
    ):
        free_free_emission(r_packet, time_explosion, numba_plasma)

    elif (
        montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED
        and transition_type == MacroAtomTransitionType.BF_EMISSION
    ):
        bound_free_emission(
            r_packet, time_explosion, numba_plasma, transition_id
        )
    elif (
        montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED
        and transition_type == MacroAtomTransitionType.BF_COOLING
    ):
        bf_cooling(r_packet, time_explosion, numba_plasma)

    elif (
        montecarlo_configuration.CONTINUUM_PROCESSES_ENABLED
        and transition_type == MacroAtomTransitionType.ADIABATIC_COOLING
    ):
        adiabatic_cooling(r_packet)

    elif transition_type == MacroAtomTransitionType.BB_EMISSION:
        line_emission(r_packet, transition_id, time_explosion, numba_plasma)
    else:
        raise Exception("No Interaction Found!")


@njit(**njit_dict_no_parallel)
def bf_cooling(r_packet, time_explosion, numba_plasma):
    """
    Bound-Free Cooling - Determine and run bf emission from cooling

    Parameters
    ----------
    r_packet : tardis.montecarlo.montecarlo_numba.r_packet.RPacket
    time_explosion : float
    numba_plasma : tardis.montecarlo.montecarlo_numba.numba_interface.NumbaPlasma
    """

    fb_cooling_prob = numba_plasma.p_fb_deactivation[
        :, r_packet.current_shell_id
    ]
    p = fb_cooling_prob[0]
    i = 0
    zrand = np.random.random()
    while p <= zrand:  # Can't search-sorted this because it's not cumulative
        i += 1
        p += fb_cooling_prob[i]
    continuum_idx = i
    bound_free_emission(r_packet, time_explosion, numba_plasma, continuum_idx)


@njit(**njit_dict_no_parallel)
def adiabatic_cooling(r_packet):
    """
    Adiabatic cooling - equivalent to destruction of the packet

    Parameters
    ----------
    r_packet: tardis.montecarlo.montecarlo_numba.r_packet.RPacket
    """

    r_packet.status = PacketStatus.ADIABATIC_COOLING


@njit(**njit_dict_no_parallel)
def get_current_line_id(nu, line_list):
    """
    Get the line id corresponding to a frequency nu in a line list

    Parameters
    ----------
    nu : float
    line_list : np.ndarray
    """

    # Note: Since this reverses the array,
    # it may be faster to just write our own reverse-binary-search

    reverse_line_list = line_list[::-1]
    number_of_lines = len(line_list)
    line_id = number_of_lines - np.searchsorted(reverse_line_list, nu)
    return line_id


@njit(**njit_dict_no_parallel)
def free_free_emission(r_packet, time_explosion, numba_plasma):
    """
    Free-Free emission - set the frequency from electron-ion interaction

    Parameters
    ----------
    r_packet : tardis.montecarlo.montecarlo_numba.r_packet.RPacket
    time_explosion : float
    numba_plasma : tardis.montecarlo.montecarlo_numba.numba_interface.NumbaPlasma
    """

    inverse_doppler_factor = get_inverse_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion
    )
    comov_nu = sample_nu_free_free(numba_plasma, r_packet.current_shell_id)
    r_packet.nu = comov_nu * inverse_doppler_factor
    current_line_id = get_current_line_id(comov_nu, numba_plasma.line_list_nu)
    r_packet.next_line_id = current_line_id

    if montecarlo_configuration.full_relativity:
        r_packet.mu = angle_aberration_CMF_to_LF(
            r_packet, time_explosion, r_packet.mu
        )


@njit(**njit_dict_no_parallel)
def bound_free_emission(r_packet, time_explosion, numba_plasma, continuum_id):
    """
    Bound-Free emission - set the frequency from photo-ionization

    Parameters
    ----------
    r_packet : tardis.montecarlo.montecarlo_numba.r_packet.RPacket
    time_explosion : float
    numba_plasma : tardis.montecarlo.montecarlo_numba.numba_interface.NumbaPlasma
    continuum_id : int
    """

    inverse_doppler_factor = get_inverse_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion
    )

    comov_nu = sample_nu_free_bound(
        numba_plasma, r_packet.current_shell_id, continuum_id
    )
    r_packet.nu = comov_nu * inverse_doppler_factor
    current_line_id = get_current_line_id(comov_nu, numba_plasma.line_list_nu)
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
    temp_doppler_factor = get_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion
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
        comov_nu = r_packet.nu * old_doppler_factor  # Is this necessary?
        r_packet.nu = comov_nu * inverse_new_doppler_factor
        activation_level_id = numba_plasma.line2macro_level_upper[
            r_packet.next_line_id
        ]
        macro_atom_event(
            activation_level_id, r_packet, time_explosion, numba_plasma
        )


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
