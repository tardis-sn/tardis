from copy import deepcopy

import pytest
import numpy as np

from tardis.simulation import Simulation
from tardis.montecarlo.montecarlo_numba import RPacket, PacketCollection


from tardis.montecarlo.montecarlo_numba.numba_interface import (
    numba_plasma_initialize, NumbaModel, Estimators, VPacketCollection)

#from tardis.montecarlo.montecarlo_numba.r_packet import (trace_packet, get_doppler_factor, get_doppler_factor_partial_relativity)

#from tardis.montecarlo import montecarlo_configuration

@pytest.fixture(scope='package')
def nb_simulation_verysimple(config_verysimple, atomic_dataset):
    atomic_data = deepcopy(atomic_dataset)
    sim = Simulation.from_config(config_verysimple, atom_data=atomic_data)
    sim.iterate(10)
    return sim


@pytest.fixture()
def verysimple_collection(nb_simulation_verysimple):
    runner = nb_simulation_verysimple.runner
    return PacketCollection(
        runner.input_nu, runner.input_mu, runner.input_energy,
        runner._output_nu, runner._output_energy)


@pytest.fixture(scope='package')
def verysimple_numba_plasma(nb_simulation_verysimple):
    return numba_plasma_initialize(nb_simulation_verysimple.plasma,
                                   line_interaction_type='macroatom')


@pytest.fixture(scope='package')
def verysimple_numba_model(nb_simulation_verysimple):
    runner = nb_simulation_verysimple.runner
    model = nb_simulation_verysimple.model
    return NumbaModel(runner.r_inner_cgs, runner.r_outer_cgs,
                             model.time_explosion.to('s').value)


@pytest.fixture(scope='package')
def verysimple_estimators(nb_simulation_verysimple):
    runner = nb_simulation_verysimple.runner
    return Estimators(runner.j_estimator, runner.nu_bar_estimator,
                      runner.j_blue_estimator, runner.Edotlu_estimator)


@pytest.fixture(scope='package')
def verysimple_vpacket_collection(nb_simulation_verysimple):
    spectrum_frequency = nb_simulation_verysimple.runner.spectrum_frequency.value
    return VPacketCollection(spectrum_frequency=spectrum_frequency,
                             number_of_vpackets=0,
                             v_packet_spawn_start_frequency=0,
                             v_packet_spawn_end_frequency=np.inf,
                             temporary_v_packet_bins=20000)


@pytest.fixture(scope='package')
def verysimple_packet_collection(nb_simulation_verysimple):
    runner = nb_simulation_verysimple.runner
    return PacketCollection(runner.input_nu, runner.input_mu,
                                         runner.input_energy,
                                         runner._output_nu,
                                         runner._output_energy)