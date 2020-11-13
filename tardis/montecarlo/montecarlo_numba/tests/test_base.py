import pytest
import pandas as pd
import os
import numpy.testing as npt
import numpy as np
from copy import deepcopy
from astropy import units as u

from tardis.montecarlo import (
    montecarlo_configuration as montecarlo_configuration,
)
import tardis.montecarlo.montecarlo_numba.base as base
from tardis.simulation import Simulation
import tardis.montecarlo.montecarlo_numba.r_packet as r_packet
import tardis.montecarlo.montecarlo_numba.single_packet_loop as spl
from tardis.montecarlo.montecarlo_numba.numba_interface import (
    PacketCollection,
    VPacketCollection,
    NumbaModel,
    numba_plasma_initialize,
    Estimators,
    configuration_initialize,
)


@pytest.mark.xfail(reason="To be implemented")
def test_montecarlo_radial1d():
    assert False


def test_montecarlo_main_loop(
    config_verysimple,
    atomic_dataset,
    tardis_ref_path,
    tmpdir,
    set_seed_fixture,
    random_call_fixture,
):

    montecarlo_configuration.LEGACY_MODE_ENABLED = True

    # Load C data from refdata
    C_fname = os.path.join(tardis_ref_path, "montecarlo_1e5_compare_data.h5")
    expected_nu = pd.read_hdf(
        C_fname, key="/simulation/runner/output_nu"
    ).values
    expected_energy = pd.read_hdf(
        C_fname, key="/simulation/runner/output_energy"
    ).values
    expected_nu_bar_estimator = pd.read_hdf(
        C_fname, key="/simulation/runner/nu_bar_estimator"
    ).values
    expected_j_estimator = pd.read_hdf(
        C_fname, key="/simulation/runner/j_estimator"
    ).values

    # Setup model config from verysimple
    atomic_data = deepcopy(atomic_dataset)
    config_verysimple.montecarlo.last_no_of_packets = 1e5
    config_verysimple.montecarlo.no_of_virtual_packets = 0
    config_verysimple.montecarlo.iterations = 1
    config_verysimple.montecarlo.single_packet_seed = 0
    del config_verysimple["config_dirname"]

    sim = Simulation.from_config(config_verysimple, atom_data=atomic_data)

    # Init model
    numba_plasma = numba_plasma_initialize(
        sim.plasma, line_interaction_type="macroatom"
    )

    runner = sim.runner
    model = sim.model

    runner._initialize_geometry_arrays(model)
    runner._initialize_estimator_arrays(numba_plasma.tau_sobolev.shape)
    runner._initialize_packets(model.t_inner.value, 100000, 0)

    # Init parameters
    montecarlo_configuration.v_packet_spawn_start_frequency = (
        runner.virtual_spectrum_spawn_range.end.to(
            u.Hz, equivalencies=u.spectral()
        ).value
    )
    montecarlo_configuration.v_packet_spawn_end_frequency = (
        runner.virtual_spectrum_spawn_range.start.to(
            u.Hz, equivalencies=u.spectral()
        ).value
    )
    montecarlo_configuration.temporary_v_packet_bins = 20000
    montecarlo_configuration.full_relativity = runner.enable_full_relativity
    montecarlo_configuration.single_packet_seed = 0

    # Init packet collection from runner
    packet_collection = PacketCollection(
        runner.input_nu,
        runner.input_mu,
        runner.input_energy,
        runner._output_nu,
        runner._output_energy,
    )

    # Init model from runner
    numba_model = NumbaModel(
        runner.r_inner_cgs,
        runner.r_outer_cgs,
        model.time_explosion.to("s").value,
    )

    # Init estimators from runner
    estimators = Estimators(
        runner.j_estimator,
        runner.nu_bar_estimator,
        runner.j_blue_estimator,
        runner.Edotlu_estimator,
    )

    # Empty vpacket collection
    vpacket_collection = VPacketCollection(
        0, np.array([0, 0], dtype=np.float64), 0, np.inf, 0, 0
    )

    # output arrays
    output_nus = np.empty_like(packet_collection.packets_output_nu)
    output_energies = np.empty_like(packet_collection.packets_output_nu)

    # IMPORTANT: seeds RNG state within JIT
    seed = 23111963
    set_seed_fixture(seed)
    for i in range(len(packet_collection.packets_input_nu)):
        # Generate packet
        packet = r_packet.RPacket(
            numba_model.r_inner[0],
            packet_collection.packets_input_mu[i],
            packet_collection.packets_input_nu[i],
            packet_collection.packets_input_energy[i],
            seed,
            i,
            0,
        )

        # Loop packet
        spl.single_packet_loop(
            packet, numba_model, numba_plasma, estimators, vpacket_collection
        )
        output_nus[i] = packet.nu
        if packet.status == r_packet.PacketStatus.REABSORBED:
            output_energies[i] = -packet.energy
        elif packet.status == r_packet.PacketStatus.EMITTED:
            output_energies[i] = packet.energy

        # RNG to match C
        random_call_fixture()

    packet_collection.packets_output_energy[:] = output_energies[:]
    packet_collection.packets_output_nu[:] = output_nus[:]

    actual_energy = packet_collection.packets_output_energy
    actual_nu = packet_collection.packets_output_nu
    actual_nu_bar_estimator = estimators.nu_bar_estimator
    actual_j_estimator = estimators.j_estimator

    # Compare
    npt.assert_allclose(
        actual_nu_bar_estimator, expected_nu_bar_estimator, rtol=1e-13
    )
    npt.assert_allclose(actual_j_estimator, expected_j_estimator, rtol=1e-13)
    npt.assert_allclose(actual_energy, expected_energy, rtol=1e-13)
    npt.assert_allclose(actual_nu, expected_nu, rtol=1e-13)
