from astropy import constants as const
import numpy as np
import pandas as pd


P_INTERNAL_UP = 1
P_INTERNAL_DOWN = 0
P_EMISSION_DOWN = -1


def line_transition_internal_up(
    line_f_lus,
    line_nus,
    energies_lower,
    mean_intensities_blue_wing,
    beta_sobolevs,
    stimulated_emission_factors,
    transition_a_i_l_u_array,
    line_ids,
):
    """
    Calculate internal upward transition probabilities for line transitions in macro atoms.

    This function computes the probability of internal upward transitions between
    energy levels in macro atoms due to line transitions. It calculates the
    transition probabilities and creates metadata for tracking the transitions.

    Parameters
    ----------
    line_f_lus : array_like
        Oscillator strengths for line transitions from lower to upper levels.
    line_nus : array_like
        Frequencies of the line transitions in Hz.
    energies_lower : array_like
        Energy values of the lower levels in the transitions.
    mean_intensities_blue_wing : array_like
        Mean radiation field intensities at the blue wing of the lines.
    beta_sobolevs : array_like
        Sobolev escape probabilities for the line transitions.
    stimulated_emission_factors : array_like
        Factors accounting for stimulated emission in the transitions.
    transition_a_i_l_u_array : array_like
        2D array containing atomic number, ion number, lower level, and upper level
        indices for each transition. Shape should be (n_transitions, 4).
    line_ids : array_like
        Unique identifiers for each line transition.

    Returns
    -------
    p_internal_up : pandas.DataFrame
        DataFrame containing the calculated internal upward transition probabilities
        with source information, indexed by transition.
    internal_up_metadata : pandas.DataFrame
        DataFrame containing metadata for the transitions including transition line IDs,
        source and destination level information, transition type, and line indices.

    Notes
    -----
    The transition probability is calculated using the formula:
    P = (f_lu / (h * nu)) * stimulated_emission_factor * mean_intensity * beta * E_lower

    The function uses the constant P_INTERNAL_UP to mark the transition type.
    """

    p_internal_up = (
        line_f_lus
        / (const.h.cgs.value * line_nus)
        * stimulated_emission_factors
        * mean_intensities_blue_wing
        * beta_sobolevs
        * energies_lower
    )
    p_internal_up["source"] = [
        tuple(col) for col in transition_a_i_l_u_array[:, [0, 1, 2]]
    ]

    p_internal_up["destination"] = [
        tuple(col) for col in transition_a_i_l_u_array[:, [0, 1, 3]]
    ]
    p_internal_up["transition_type"] = P_INTERNAL_UP

    p_internal_up["transition_line_id"] = line_ids

    p_internal_up["transition_line_idx"] = range(len(line_ids))

    internal_up_metadata = p_internal_up[
        [
            "transition_line_id",
            "source",
            "destination",
            "transition_type",
            "transition_line_idx",
        ]
    ]

    p_internal_up = p_internal_up.drop(
        columns=[
            "destination",
            "transition_type",
            "transition_line_id",
            "transition_line_idx",
        ]
    )

    return p_internal_up, internal_up_metadata


def line_transition_internal_down(
    line_f_uls,
    line_nus,
    energies_lower,
    beta_sobolevs,
    transition_a_i_l_u_array,
    line_ids,
):
    """
    Calculate internal downward transition probabilities for line transitions.

    This function computes the probability of internal downward transitions
    in macro atoms for line transitions, based on oscillator strengths,
    frequencies, and other atomic parameters.

    Parameters
    ----------
    line_f_uls : array_like
        Oscillator strengths for line transitions.
    line_nus : array_like
        Frequencies of line transitions.
    energies_lower : array_like
        Energies of the lower levels involved in transitions.
    beta_sobolevs : array_like
        Sobolev beta factors for the transitions.
    transition_a_i_l_u_array : array_like
        Array containing atomic number, ion number, lower level, and upper level
        indices for each transition. Shape should be (n_transitions, 4).
    line_ids : array_like
        Identifiers for each line transition.

    Returns
    -------
    p_internal_down : pandas.DataFrame
        DataFrame containing the calculated internal downward transition
        probabilities with source information.
    internal_down_metadata : pandas.DataFrame
        DataFrame containing metadata for the transitions including
        transition line IDs, source and destination level information,
        transition type, and transition line indices.

    Notes
    -----
    The internal downward transition probability is calculated using the formula:
    p = 2 * nu^2 * f_ul / c^2 * beta * E_lower

    The function uses the constant P_INTERNAL_DOWN to mark the transition type.
    """
    p_internal_down = (
        2
        * line_nus**2
        * line_f_uls
        / const.c.cgs.value**2
        * beta_sobolevs
        * energies_lower
    )
    p_internal_down["source"] = [
        tuple(col) for col in transition_a_i_l_u_array[:, [0, 1, 3]]
    ]

    p_internal_down["destination"] = [
        tuple(col) for col in transition_a_i_l_u_array[:, [0, 1, 2]]
    ]
    p_internal_down["transition_type"] = P_INTERNAL_DOWN

    p_internal_down["transition_line_id"] = line_ids

    p_internal_down["transition_line_idx"] = range(len(line_ids))

    internal_down_metadata = p_internal_down[
        [
            "transition_line_id",
            "source",
            "destination",
            "transition_type",
            "transition_line_idx",
        ]
    ]

    p_internal_down = p_internal_down.drop(
        columns=[
            "destination",
            "transition_type",
            "transition_line_id",
            "transition_line_idx",
        ]
    )
    return p_internal_down, internal_down_metadata


def line_transition_emission_down(
    line_f_uls,
    line_nus,
    energies_upper,
    energies_lower,
    beta_sobolevs,
    transition_a_i_l_u_array,
    line_ids,
):
    """
    Calculate emission down transition probabilities for line transitions.

    This function computes the probability of emission down transitions for
    atomic line transitions using oscillator strengths, frequencies, energy
    differences, and Sobolev optical depths.

    Parameters
    ----------
    line_f_uls : array_like
        Oscillator strengths for the line transitions.
    line_nus : array_like
        Frequencies of the line transitions in Hz.
    energies_upper : array_like
        Energy values of the upper levels in the transitions.
    energies_lower : array_like
        Energy values of the lower levels in the transitions.
    beta_sobolevs : array_like
        Sobolev escape probabilities for the transitions.
    transition_a_i_l_u_array : array_like
        Array containing atomic number, ion number, lower level, and upper level
        indices for each transition. Shape should be (n_transitions, 4).
    line_ids : array_like
        Unique identifiers for each line transition.

    Returns
    -------
    p_emission_down : pandas.DataFrame
        DataFrame containing the calculated emission down probabilities with
        source information. Contains columns for the probability values and
        'source' tuples.
    emission_down_metadata : pandas.DataFrame
        DataFrame containing metadata for the transitions including transition
        line IDs, source and destination level information, transition type,
        and transition line indices.

    Notes
    -----
    The emission down probability is calculated using the formula:
    P = 2 * nu^2 * f_ul / c^2 * beta * (E_upper - E_lower)

    The function uses the constant P_EMISSION_DOWN to mark the transition type.
    """
    p_emission_down = (
        2
        * line_nus**2
        * line_f_uls
        / const.c.cgs.value**2
        * beta_sobolevs
        * (energies_upper - energies_lower)
    )
    p_emission_down["source"] = [
        tuple(col) for col in transition_a_i_l_u_array[:, [0, 1, 3]]
    ]

    p_emission_down["destination"] = [
        tuple(col) for col in transition_a_i_l_u_array[:, [0, 1, 2]]
    ]

    p_emission_down["transition_type"] = P_EMISSION_DOWN

    p_emission_down["transition_line_id"] = line_ids

    p_emission_down["transition_line_idx"] = range(len(line_ids))

    emission_down_metadata = p_emission_down[
        [
            "transition_line_id",
            "source",
            "destination",
            "transition_type",
            "transition_line_idx",
        ]
    ]

    p_emission_down = p_emission_down.drop(
        columns=[
            "destination",
            "transition_type",
            "transition_line_id",
            "transition_line_idx",
        ]
    )
    return p_emission_down, emission_down_metadata



def set_index(p, photo_ion_idx, transition_type=0, reverse=True):
        idx = photo_ion_idx.loc[p.index]
        transition_type = transition_type * np.ones_like(
            idx.destination_level_idx
        )
        transition_type = pd.Series(transition_type, name="transition_type")
        idx_arrays = [idx.source_level_idx, idx.destination_level_idx]
        if reverse:
            idx_arrays = idx_arrays[::-1]
        idx_arrays.append(transition_type)
        index = pd.MultiIndex.from_arrays(idx_arrays)
        if reverse:
            index.names = index.names[:-1][::-1] + [index.names[-1]]
        p = p.set_index(index, drop=True)
        return p

def continuum_transition_recombination(alpha_sp, nu_i, energy_i, photoionization_index):
    '''Unnormalized probabilities of radiative recombination
    alpha_sp : pandas.Series, dtype float
        Rate coefficient for spontaneous recombination from `k` to level `i`
    nu_i : pandas.Series, dtype float
        Threshold frequencies for ionization
    energy_i : pandas.Series, dtype float
        Energies of levels with bound-free transitions. Needed to calculate
        for example internal transition probabilities in the macro atom scheme.
    All p_internals (u/l) go as R_i(u/l) * e_i / D_i
    All p_deactivations go as R * (delta_e)/D_i'''

    p_recomb_deactivation = alpha_sp.multiply(nu_i, axis=0) * const.h.cgs.value # nu_i is ionization threshold, so I guess delta e is just nu_i * h?
    p_recomb_deactivation = set_index(
        p_recomb_deactivation, photoionization_index, transition_type=P_EMISSION_DOWN
    )

    p_recomb_internal = alpha_sp.multiply(energy_i, axis=0)
    p_recomb_internal = set_index(
        p_recomb_internal, photoionization_index, transition_type=P_INTERNAL_DOWN
    )
    p_recombination = pd.concat([p_recomb_deactivation, p_recomb_internal])

    recombination_metadata = p_recombination[
        [
            "transition_line_id",
            "source",
            "destination",
            "transition_type",
            "transition_line_idx",
        ]
    ]

    p_recombination = p_recombination.drop(
        columns=[
            "destination",
            "transition_type",
            "transition_line_id",
            "transition_line_idx",
        ]
    )

    return p_recombination, recombination_metadata


def continuum_transition_photoionization(gamma_corr, energy_i, photo_ion_idx):
    '''
    Photoionization probability unnormalized
        gamma_corr : pandas.Series, dtype float
            Corrected photoionization rate coefficeint from level `i` to `k`
    This method may need an additional transition type
     '''
    p_photoionization = gamma_corr.multiply(energy_i, axis=0)
    p_photoionization = set_index(p_photoionization, photo_ion_idx, reverse=False, transition_type=P_INTERNAL_DOWN)

    photoionization_metadata = p_photoionization[
        [
            "transition_line_id",
            "source",
            "destination",
            "transition_type",
            "transition_line_idx",
        ]
    ]

    p_photoionization = p_photoionization.drop(
        columns=[
            "destination",
            "transition_type",
            "transition_line_id",
            "transition_line_idx",
        ]
    )

    return p_photoionization, photoionization_metadata




def continuum_transition_collisional(
    self,
    coll_exc_coeff,
    coll_deexc_coeff,
    yg_idx,
    electron_densities,
    delta_E_yg,
    atomic_data,
    level_number_density,
):
    p_deexc_deactivation = (coll_deexc_coeff * electron_densities).multiply(
        delta_E_yg.values, axis=0
    )
    p_deexc_deactivation = self.set_index(
        p_deexc_deactivation, yg_idx, reverse=True
    )
    p_deexc_deactivation = p_deexc_deactivation.groupby(level=[0]).sum()
    index_dd = pd.MultiIndex.from_product(
        [p_deexc_deactivation.index.values, ["k"], [0]],
        names=list(yg_idx.columns) + ["transition_type"],
    )
    p_deexc_deactivation = p_deexc_deactivation.set_index(index_dd)

    level_lower_index = coll_deexc_coeff.index.droplevel(
        "level_number_upper"
    )
    energy_lower = atomic_data.levels.energy.loc[level_lower_index]
    p_deexc_internal = (coll_deexc_coeff * electron_densities).multiply(
        energy_lower.values, axis=0
    )
    p_deexc_internal = self.set_index(
        p_deexc_internal, yg_idx, transition_type=0, reverse=True
    )

    p_exc_internal = (coll_exc_coeff * electron_densities).multiply(
        energy_lower.values, axis=0
    )
    p_exc_internal = self.set_index(
        p_exc_internal, yg_idx, transition_type=0, reverse=False
    )
    p_exc_cool = (coll_exc_coeff * electron_densities).multiply(
        delta_E_yg.values, axis=0
    )
    p_exc_cool = (
        p_exc_cool * level_number_density.loc[level_lower_index].values
    )
    p_exc_cool = self.set_index(p_exc_cool, yg_idx, reverse=False)
    p_exc_cool = p_exc_cool.groupby(level="destination_level_idx").sum()
    exc_cool_index = pd.MultiIndex.from_product(
        [["k"], p_exc_cool.index.values, [0]],
        names=list(yg_idx.columns) + ["transition_type"],
    )
    p_exc_cool = p_exc_cool.set_index(exc_cool_index)
    p_coll = pd.concat(
        [p_deexc_deactivation, p_deexc_internal, p_exc_internal, p_exc_cool]
    )
    return p_coll


