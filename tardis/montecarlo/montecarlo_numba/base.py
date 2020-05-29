from numba import prange, njit, jit
import logging
import numpy as np

# from tardis.montecarlo import __init__
from tardis.montecarlo.montecarlo_numba.r_packet import RPacket, PacketStatus
from tardis.montecarlo.montecarlo_numba.numba_interface import (
    PacketCollection, VPacketCollection, NumbaModel, numba_plasma_initialize,
    Estimators, MonteCarloConfiguration, configuration_initialize)

from tardis.montecarlo.montecarlo_numba.single_packet_loop import (
    single_packet_loop)
from tardis.montecarlo.montecarlo_numba import njit_dict

logger = logging.getLogger('tardis-montecarlo-numba-debug-logger')


def montecarlo_radial1d(model, plasma, runner, montecarlo_configuration):
    packet_collection = PacketCollection(
        runner.input_nu, runner.input_mu, runner.input_energy,
        runner._output_nu, runner._output_energy
    )

    numba_model = NumbaModel(runner.r_inner_cgs, runner.r_outer_cgs,
                             model.time_explosion.to('s').value)
    numba_plasma = numba_plasma_initialize(plasma)
    estimators = Estimators(runner.j_estimator, runner.nu_bar_estimator,
                            runner.j_b_lu_estimator, runner.edot_lu_estimator)

    v_packets_energy_hist = montecarlo_main_loop(
        packet_collection, numba_model, numba_plasma, estimators,
        runner.spectrum_frequency.value, montecarlo_configuration)
    
    runner._montecarlo_virtual_luminosity.value[:] = v_packets_energy_hist


@jit(**njit_dict, nogil=True)
def montecarlo_main_loop(packet_collection, numba_model, numba_plasma,
                         estimators, spectrum_frequency,
                         montecarlo_configuration):
    """
    This is the main loop of the MonteCarlo routine that generates packets 
    and sends them through the ejecta. 

    Parameters
    ----------
    storage_model : [type]
        [description]
    """

    output_nus = np.empty_like(packet_collection.packets_output_nu)
    output_energies = np.empty_like(packet_collection.packets_output_nu)

    v_packets_energy_hist = np.zeros_like(spectrum_frequency)
    delta_nu = spectrum_frequency[1] - spectrum_frequency[0]

    for i in prange(len(output_nus)):
        r_packet = RPacket(numba_model.r_inner[0],
                           packet_collection.packets_input_mu[i],
                           packet_collection.packets_input_nu[i],
                           packet_collection.packets_input_energy[i],
                           i)
        np.random.seed(r_packet.index)
        vpacket_collection = VPacketCollection(
            spectrum_frequency, montecarlo_configuration.number_of_vpackets,
            montecarlo_configuration.temporary_v_packet_bins)
        single_packet_loop(r_packet, numba_model, numba_plasma, estimators,
                           vpacket_collection,
                           montecarlo_configuration)

        output_nus[i] = r_packet.nu

        if r_packet.status == PacketStatus.REABSORBED:
            output_energies[i] = -r_packet.energy
        elif r_packet.status == PacketStatus.EMITTED:
            output_energies[i] = r_packet.energy
        
        vpackets_nu = vpacket_collection.nus[:vpacket_collection.idx]
        vpackets_energy = vpacket_collection.energies[:vpacket_collection.idx]

        v_packets_idx = np.floor((vpackets_nu - spectrum_frequency[0]) /
                                 delta_nu).astype(np.int64)
        for j, idx in enumerate(v_packets_idx):
            if ((vpackets_nu[j] < spectrum_frequency[0]) or
                    (vpackets_nu[j] > spectrum_frequency[-1])):
                continue
            v_packets_energy_hist[idx] += vpackets_energy[j]

    packet_collection.packets_output_energy[:] = output_energies[:]
    packet_collection.packets_output_nu[:] = output_nus[:]
    
    return v_packets_energy_hist

def log_decorator(func):
    """
    Decorator to log functions while in debug mode, i.e., when
    `debug_montecarlo` is True in the config. Works for
    `@jit'd and `@njit`'d functions, but with a significant speed
    penalty.

    Questions:
        - stdout or print to file?

    TODO: How do log *args?
    TODO: How to pass kwargs to @jit, @njit?
    TODO: in nopython mode: do I need a context manager?
    TODO: Buffer?
    TODO: make numpy docstring.
    TODO: have this know debug_mode from the config.

    :param func: function to be logged.
    :return: either the function itself, if debug_mode is true, or
    """
    if __init__.DEBUG_MODE:
        logger.setLevel(logging.DEBUG)
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.DEBUG)
        console_formatter = logging.Formatter(
            '%(name)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)

        def wrapper(*args, **kwargs):
            logger.debug(f'Func: {func}. Input: {(args, kwargs)}')
            result = func(*args)
            logger.debug(f'Output: {result}.')
            return result

        return wrapper
    else:
        return func