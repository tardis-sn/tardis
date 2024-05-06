import pytest
import numpy.testing as npt
import numpy as np
import tardis.montecarlo.montecarlo_numba.interaction as interaction
from tardis.montecarlo.montecarlo_numba.numba_interface import (
    LineInteractionType,
)


def test_thomson_scatter(packet, verysimple_numba_model):
    init_mu = packet.mu
    init_nu = packet.nu
    init_energy = packet.energy
    time_explosion = verysimple_numba_model.time_explosion

    interaction.thomson_scatter(packet, time_explosion, False)

    assert np.abs(packet.mu - init_mu) > 1e-7
    assert np.abs(packet.nu - init_nu) > 1e-7
    assert np.abs(packet.energy - init_energy) > 1e-7


@pytest.mark.parametrize(
    "line_interaction_type",
    [
        LineInteractionType.SCATTER,
        LineInteractionType.DOWNBRANCH,
        LineInteractionType.MACROATOM,
    ],
)
def test_line_scatter(
    line_interaction_type,
    packet,
    verysimple_numba_model,
    verysimple_opacity_state,
):
    init_mu = packet.mu
    init_nu = packet.nu
    init_energy = packet.energy
    full_relativity = False
    packet.initialize_line_id(
        verysimple_opacity_state, verysimple_numba_model, full_relativity
    )
    time_explosion = verysimple_numba_model.time_explosion

    interaction.line_scatter(
        packet,
        time_explosion,
        line_interaction_type,
        verysimple_opacity_state,
        continuum_processes_enabled=False,
        enable_full_relativity=False,
    )

    assert np.abs(packet.mu - init_mu) > 1e-7
    assert np.abs(packet.nu - init_nu) > 1e-7
    assert np.abs(packet.energy - init_energy) > 1e-7


@pytest.mark.parametrize(
    ["test_packet", "expected"],
    [
        (
            {
                "mu": 0.8599443103322428,
                "emission_line_id": 1000,
                "energy": 0.9114437898710559,
                "nu": 0.0,
            },
            {"mu": 0.8599443103322428, "energy": 0.9114437898710559},
        ),
        (
            {
                "mu": -0.6975116557422458,
                "emission_line_id": 2000,
                "energy": 0.8803098648913266,
            },
            {"mu": -0.6975116557422458, "energy": 0.8803098648913266},
        ),
        (
            {
                "mu": -0.7115661419975774,
                "emission_line_id": 0,
                "energy": 0.8800385929341252,
            },
            {"mu": -0.7115661419975774, "energy": 0.8800385929341252},
        ),
    ],
)
def test_line_emission(
    packet,
    verysimple_numba_model,
    verysimple_opacity_state,
    test_packet,
    expected,
):
    emission_line_id = test_packet["emission_line_id"]
    packet.mu = test_packet["mu"]
    packet.energy = test_packet["energy"]
    full_relativity = False
    packet.initialize_line_id(
        verysimple_opacity_state,
        verysimple_numba_model,
        full_relativity,
    )

    time_explosion = verysimple_numba_model.time_explosion

    interaction.line_emission(
        packet,
        emission_line_id,
        time_explosion,
        verysimple_opacity_state,
        full_relativity,
    )

    assert packet.next_line_id == emission_line_id + 1
    npt.assert_almost_equal(packet.mu, expected["mu"])
    npt.assert_almost_equal(packet.energy, expected["energy"])
