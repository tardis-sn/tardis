import photon
import constants
import numpy as np
import montecarlo_multizone
import time
import model
import os
import logging
import synspec
import texttable



# Adding logging support
logger = logging.getLogger(__name__)


def run_multizone(config_dict, atomic_model):
    line2level = atomic_model.line_list['global_level_id_upper'] - 1
    atomic_model.macro_atom.read_nus(atomic_model.line_list['nu'])

    line_interaction_type = config_dict.get('line_interaction_type', 'macro')
    logger.debug('Line interaction type set to \'%s\'' % line_interaction_type)

    if config_dict['config_type'] == 'uniform_w7':
        current_model = model.MultiZoneRadial.from_lucy99(config_dict['v_inner'], config_dict['abundances'],
            config_dict['time_exp'], config_dict['no_of_shells'])

    elif config_dict['config_type'] == 'shell':
        current_model = model.MultiZoneRadial(config_dict['velocities'], config_dict['densities'],
            config_dict['abundances'], config_dict['time_exp'], config_dict['ws'])

    surface_inner = 4 * np.pi * current_model.r_inner[0] ** 2
    volume = (4. / 3) * np.pi * (current_model.r_outer ** 3 - current_model.r_inner ** 3)
    current_model.set_atomic_model(atomic_model)

    current_model.initialize_plasmas(config_dict['t_rads'])

    logger.debug('\n\n' + current_model.print_model_table())
    #TODO temperature initialized wrong. Need to calculate inner temperature not outer temperature
    t_inner = config_dict['t_rads'][0]
    if line_interaction_type == 'macro':
        do_scatter = 0

    elif line_interaction_type == 'downbranch':
        raise NotImplementedError('downbranch yet to be implemented')

    elif line_interaction_type == 'scatter':
        do_scatter = 1

    else:
        raise ValueError('Line interaction type %s not understood (allowed are macro, downbranch, scatter)'\
                         % line_interaction_type)

    i = 0
    track_ws = []
    track_t_rads = []
    track_t_inner = []

    while True:
        start_time = time.time()
        track_t_rads.append(current_model.t_rads.copy())
        track_t_inner.append(t_inner)
        track_ws.append(current_model.ws.copy())
        i += 1
        if i > config_dict['iterations']: break
        if os.path.exists('stop_file') or i == config_dict['iterations']:
            no_of_packets = config_dict['no_of_spectrum_packets']
            energy_of_packet = 1. / no_of_packets
        else:
            no_of_packets = config_dict['no_of_calibration_packets']
            energy_of_packet = 1. / no_of_packets

        #return sn_plasma



        tau_sobolevs = current_model.calculate_tau_sobolevs()

        if config_dict['exclude_ions'] is not None:
            ion_filter = np.zeros_like(atomic_model.line_list['ion']).astype(bool)
            for ion in config_dict['exclude_ions']:
                ion_filter = ion_filter | atomic_model.line_list['ion'] != ion - 1

            tau_sobolevs[0][ion_filter] = 0.

        transition_probabilities = current_model.calculate_transition_probabilities(tau_sobolevs)
        nu_input = np.sort(
            photon.random_blackbody_nu(t_inner, nu_range=(1e8 * constants.c / 1, 1e8 * constants.c / 100000.),
                size=no_of_packets))[::-1]
        mu_input = np.sqrt(np.random.random(no_of_packets))

        logger.info("Start TARDIS MonteCarlo")

        nu, energy, nu_reabsorbed, energy_reabsorbed, j_estimators, nubar_estimators =\
        montecarlo_multizone.run_simple_oned(nu_input,
            mu_input,
            atomic_model.line_list['nu'],
            tau_sobolevs,
            current_model.r_inner,
            current_model.r_outer,
            current_model.v_inner,
            current_model.electron_densities,
            energy_of_packet,
            transition_probabilities,
            current_model.atomic_model.macro_atom.transition_type_total,
            current_model.atomic_model.macro_atom.target_level_total,
            current_model.atomic_model.macro_atom.target_line_total,
            current_model.atomic_model.macro_atom.level_references,
            line2level,
            log_packets=0,
            do_scatter=do_scatter)

        new_t_rads = constants.trad_estimator_constant * nubar_estimators / j_estimators
        new_ws = constants.w_estimator_constant * j_estimators / (
            config_dict['time_of_simulation'] * (new_t_rads ** 4) * volume)

        emitted_energy_fraction = np.sum(energy[nu != 0]) / 1.
        reabsorbed_energy_fraction = np.sum(energy_reabsorbed[nu_reabsorbed != 0]) / 1.
        #TODO find out where the energy went. reabsorbed + emitted = 0.98
        #print "testing energy fraction emitted energy_fraction %s reabsorbed energy fraction %2 " % (
        #emitted_energy_fraction, reabsorbed_energy_fraction)

        new_t_inner = (config_dict['luminosity_outer'] / (
            emitted_energy_fraction * constants.sigma_sb * surface_inner )) ** .25

        if logger.getEffectiveLevel() == logging.DEBUG:
            temp_table = texttable.Texttable()
            header = ('t_rad', 'new_t_rad', 'w', 'new_w')
            temp_table.add_row(header)
            temp_table.set_deco(temp_table.HEADER | temp_table.VLINES)
            logger.debug("t_inner = %.2f new_tinner = %.2f", t_inner, new_t_inner)
            for new_t_rad, new_w, old_t_rad, old_w in zip(new_t_rads, new_ws, current_model.t_rads, current_model.ws):
                temp_table.add_row(('%.2f' % old_t_rad, '%.2f' % new_t_rad, '%.4f' % old_w, '%.4f' % new_w))
            logger.debug('\n\n' + temp_table.draw())

        current_model.ws = (new_ws + current_model.ws) * .5
        current_model.t_rads = (new_t_rads + current_model.t_rads) * .5
        t_inner = (new_t_inner + t_inner) * .5
        current_model.update_model(current_model.t_rads, current_model.ws)

        logger.info("Last iteration took %.2f s", (time.time() - start_time))
    return synspec.tardis_result(nu_input,
        energy_of_packet,
        nu,
        energy,
        nu_reabsorbed,
        energy_reabsorbed,
        track_t_rads,
        track_ws,
        track_t_inner)
