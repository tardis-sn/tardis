import numpy as np
import numpy.testing as npt
import pytest

from tardis.transport.montecarlo import RPacket
from tardis.transport.montecarlo.configuration.base import (
    MonteCarloConfiguration,
)
from tardis.transport.montecarlo.estimators.estimators_bulk import (
    init_estimators_bulk,
)
from tardis.transport.montecarlo.estimators.estimators_line import (
    init_estimators_line,
)
from tardis.transport.montecarlo.interaction_events import LineInteractionType
from tardis.transport.montecarlo.modes.classic.packet_propagation import (
    packet_propagation,
)
from tardis.transport.montecarlo.packets.trackers.tracker_last_interaction import (
    TrackerLastInteraction,
)


def test_verysimple_single_packet_loop(
    verysimple_numba_radial_1d_geometry,
    verysimple_time_explosion,
    verysimple_opacity_state,
    verysimple_vpacket_collection,
    verysimple_packet_collection,
):
    numba_radial_1d_geometry = verysimple_numba_radial_1d_geometry
    packet_collection = verysimple_packet_collection
    vpacket_collection = verysimple_vpacket_collection
    time_explosion = verysimple_time_explosion
    opacity_state = verysimple_opacity_state

    n_cells = len(numba_radial_1d_geometry.r_inner)
    n_lines_by_n_cells_tuple = opacity_state.tau_sobolev.shape
    estimators_bulk = init_estimators_bulk(n_cells)
    estimators_line = init_estimators_line(n_lines_by_n_cells_tuple)

    rpacket_tracker = TrackerLastInteraction()
    montecarlo_configuration = MonteCarloConfiguration()
    montecarlo_configuration.LINE_INTERACTION_TYPE = LineInteractionType.MACROATOM

    i = 0
    np.random.seed(packet_collection.packet_seeds[i])
    r_packet = RPacket(
        r=numba_radial_1d_geometry.r_inner[0],
        mu=packet_collection.initial_mus[i],
        nu=packet_collection.initial_nus[i],
        energy=packet_collection.initial_energies[i],
        seed=packet_collection.packet_seeds[i],
        index=i,
    )

    packet_propagation(
        r_packet,
        numba_radial_1d_geometry,
        time_explosion,
        opacity_state,
        estimators_bulk,
        estimators_line,
        vpacket_collection,
        rpacket_tracker,
        montecarlo_configuration,
    )

    assert r_packet.status != 0
    npt.assert_almost_equal(r_packet.nu, 1405610115898994.5)
    npt.assert_almost_equal(r_packet.mu, 0.9611146425440567)
    npt.assert_almost_equal(r_packet.energy, 0.10327717505563379)


@pytest.mark.xfail(reason="To be implemented")
def test_set_packet_props_partial_relativity():
    raise AssertionError()


@pytest.mark.xfail(reason="To be implemented")
def test_set_packet_props_full_relativity():
    raise AssertionError()
