@njit(**njit_dict_no_parallel)
def trace_packet_continuum(
        r_packet,
        numba_model,
        numba_plasma,
        estimators,
        continuum
):
    """
    Traces the RPacket through the ejecta and stops when an interaction happens (heart of the calculation)

    Parameters
    ----------
    r_packet : tardis.montecarlo.montecarlo_numba.r_packet.RPacket
    numba_model : tardis.montecarlo.montecarlo_numba.numba_interface.NumbaModel
    numba_plasma : tardis.montecarlo.montecarlo_numba.numba_interface.NumbaPlasma
    estimators : tardis.montecarlo.montecarlo_numba.numba_interface.Estimators

    Returns
    -------
    """

    r_inner = numba_model.r_inner[r_packet.current_shell_id]
    r_outer = numba_model.r_outer[r_packet.current_shell_id]

    (
        distance_boundary,
        delta_shell,
    ) = calculate_distance_boundary(r_packet.r, r_packet.mu, r_inner, r_outer)

    # defining start for line interaction
    start_line_id = r_packet.next_line_id

    # defining taus
    tau_event = -np.log(np.random.random())
    tau_trace_line_combined = 0.0

    # e scattering initialization

    cur_electron_density = numba_plasma.electron_density[
        r_packet.current_shell_id
    ]
    chi_e = cur_electron_density * SIGMA_THOMSON

    # Calculating doppler factor
    doppler_factor = get_doppler_factor(
        r_packet.r, r_packet.mu, numba_model.time_explosion
    )
    comov_nu = r_packet.nu * doppler_factor
    continuum.calculate(comov_nu, r_packet.current_shell_id)
    (
        chi_bf,
        chi_bf_contributions,
        current_continua,
        x_sect_bfs,
        chi_ff,
    ) = (
            continuum.chi_bf_tot, 
            continuum.chi_bf_contributions, 
            continuum.current_continua, 
            continuum.x_sect_bfs, 
            continuum.chi_ff
    )

    chi_continuum = chi_e + chi_bf + chi_ff
    distance_continuum = tau_event / chi_continuum

    cur_line_id = start_line_id  # initializing varibale for Numba
    # - do not remove
    last_line_id = len(numba_plasma.line_list_nu) - 1
    for cur_line_id in range(start_line_id, len(numba_plasma.line_list_nu)):

        # Going through the lines
        nu_line = numba_plasma.line_list_nu[cur_line_id]
        nu_line_last_interaction = numba_plasma.line_list_nu[cur_line_id - 1]

        # Getting the tau for the next line
        tau_trace_line = numba_plasma.tau_sobolev[
            cur_line_id, r_packet.current_shell_id
        ]

        # Adding it to the tau_trace_line_combined
        tau_trace_line_combined += tau_trace_line

        # Calculating the distance until the current photons co-moving nu
        # redshifts to the line frequency
        is_last_line = cur_line_id == last_line_id

        distance_trace = calculate_distance_line(
            r_packet,
            comov_nu,
            is_last_line,
            nu_line,
            numba_model.time_explosion,
        )

        # calculating the tau continuum of how far the trace has progressed
        tau_trace_continuum = chi_continuum * distance_trace

        # calculating the trace
        tau_trace_combined = tau_trace_line_combined + tau_trace_continuum

        distance = min(distance_trace, distance_boundary, distance_continuum)

        if distance_trace != 0:

            if distance == distance_boundary:
                interaction_type = InteractionType.BOUNDARY  # BOUNDARY
                r_packet.next_line_id = cur_line_id
                break
            elif distance == distance_continuum:
                zrand = np.random.random()
                if zrand < chi_e / chi_continuum:
                    interaction_type = InteractionType.ESCATTERING
                else:
                    interaction_type = InteractionType.CONTINUUM_PROCESS
                r_packet.next_line_id = cur_line_id
                break

        # Updating the J_b_lu and E_dot_lu
        # This means we are still looking for line interaction and have not
        # been kicked out of the path by boundary or electron interaction

        update_line_estimators(
            estimators,
            r_packet,
            cur_line_id,
            distance_trace,
            numba_model.time_explosion,
        )

        if (
            tau_trace_combined > tau_event
            and not montecarlo_configuration.disable_line_scattering
        ):
            interaction_type = InteractionType.LINE  # Line
            r_packet.last_interaction_in_nu = r_packet.nu
            r_packet.last_line_interaction_in_id = cur_line_id
            r_packet.next_line_id = cur_line_id
            distance = distance_trace
            break

        # Recalculating distance_continuum using tau_event -
        # tau_trace_line_combined
        # I don't think this needs to be updated
        # since tau_event is already the result of the integral
        # from the initial line
        distance_continuum = (tau_event - tau_trace_line_combined) / (
            chi_continuum
        )

    else:  # Executed when no break occurs in the for loop
        # We are beyond the line list now and the only next thing is to see
        # if we are interacting with the boundary or electron scattering
        if cur_line_id == (len(numba_plasma.line_list_nu) - 1):
            # Treatment for last line
            cur_line_id += 1
        if distance_continuum < distance_boundary:
            distance = distance_continuum
            zrand = np.random.random()
            if zrand < chi_e / chi_continuum:
                interaction_type = InteractionType.ESCATTERING
            else:
                interaction_type = InteractionType.CONTINUUM_PROCESS
            # #print('scattering')
        else:
            distance = distance_boundary
            interaction_type = InteractionType.BOUNDARY

    # r_packet.next_line_id = cur_line_id
    update_bound_free_estimators(
        comov_nu,
        r_packet.energy * doppler_factor,
        r_packet.current_shell_id,
        distance,
        estimators,
        numba_plasma.t_electrons[r_packet.current_shell_id],
        x_sect_bfs,
        current_continua,
        numba_plasma.bf_threshold_list_nu,
    )

    return distance, interaction_type, delta_shell