from astropy import constants as const

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
    line_f_lus,
    line_nus,
    energies_lower,
    beta_sobolevs,
    transition_a_i_l_u_array,
    line_ids,
):
    p_internal_down = (
        2
        * line_nus**2
        * line_f_lus
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
    line_f_lus,
    line_nus,
    energies_upper,
    energies_lower,
    beta_sobolevs,
    transition_a_i_l_u_array,
    line_ids,
):
    p_emission_down = (
        2
        * line_nus**2
        * line_f_lus
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
