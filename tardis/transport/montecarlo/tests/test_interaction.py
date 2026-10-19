from collections.abc import Callable

import numpy as np
import numpy.testing as npt
import pytest

import tardis.transport.montecarlo.interaction_events as interaction_events
from tardis import constants as const
from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.transport.montecarlo import RPacket
from tardis.transport.montecarlo.interaction_event_callers import (
    line_scatter_event,
)
from tardis.transport.montecarlo.modes.nonhomologous.interaction_event_callers import (
    line_scatter_event as nonhomologous_line_scatter_event,
)


@pytest.mark.parametrize(
    "line_interaction_type", list(interaction_events.LineInteractionType)
)
@pytest.mark.parametrize(
    "regime",
    [
        "static_nonhomologous",
        "finite_homologous_partial",
        "finite_homologous_full",
    ],
)
def test_line_redistribution_conserves_comoving_packet_energy(
    line_interaction_type: interaction_events.LineInteractionType,
    regime: str,
    parametrized_packet: RPacket,
    classic_opacity_state: OpacityStateNumba,
    set_seed_fixture: Callable[[int], None],
) -> None:
    """Verify radiative redistribution conserves comoving energy (verification).

    Claim: SCATTER, DOWNBRANCH and MACROATOM preserve packet energy in the
    comoving frame and emit the selected BB line. Regime: exactly static
    nonhomologous geometry or finite-velocity homologous partial/full, with a
    single probability-one radiative deactivation and no cooling channels.
    Verification: E*gamma*(1-beta*mu) before/after (gamma=1 in partial), and
    emitted comoving frequency from the forced line. rtol=1e-12, atol=0 covers
    floating point transforms. Detects lost energy, a wrong BB channel, or a
    missing/reversed Doppler or aberration transform.
    """
    packet = parametrized_packet
    packet.next_line_id = 0
    opacity = classic_opacity_state
    opacity.line_list_nu[:] = [4.0e14, 3.0e14]
    opacity.transition_probabilities[:] = 1.0
    opacity.macro_block_edge_index = np.array([0, 1], dtype=np.int64)
    opacity.transition_type[:] = -1  # BB emission, with probability one.
    opacity.transition_line_id[:] = 1
    time_explosion = (
        5.2e5  # s; finite, never an infinite-time static surrogate.
    )
    beta = packet.r / time_explosion / const.c.cgs.value
    full = regime == "finite_homologous_full"
    if regime == "static_nonhomologous":
        beta = 0.0
    gamma = 1.0 / np.sqrt(1.0 - beta**2) if full else 1.0
    initial_comoving_energy = packet.energy * gamma * (1.0 - beta * packet.mu)
    expected_line_id = (
        0
        if line_interaction_type
        == interaction_events.LineInteractionType.SCATTER
        else 1
    )
    set_seed_fixture(1963)
    if regime == "static_nonhomologous":
        geometry = NumbaRadial1DGeometry(
            np.array([7.0e14]), np.array([8.0e14]), np.zeros(1), np.zeros(1)
        )
        nonhomologous_line_scatter_event(
            packet, geometry, line_interaction_type, opacity, False
        )
    else:
        line_scatter_event(
            packet, time_explosion, line_interaction_type, opacity, full
        )
    final_doppler = gamma * (1.0 - beta * packet.mu)
    npt.assert_allclose(
        packet.energy * final_doppler,
        initial_comoving_energy,
        rtol=1.0e-12,
        atol=0,
    )
    npt.assert_allclose(
        packet.nu * final_doppler,
        opacity.line_list_nu[expected_line_id],
        rtol=1.0e-12,
        atol=0,
    )
    assert packet.next_line_id == expected_line_id + 1


def test_thomson_scatter(packet, verysimple_time_explosion):
    init_mu = packet.mu
    init_nu = packet.nu
    init_energy = packet.energy
    time_explosion = verysimple_time_explosion

    interaction_events.thomson_scatter(packet, time_explosion, False)

    assert np.abs(packet.mu - init_mu) > 1e-7
    assert np.abs(packet.nu - init_nu) > 1e-7
    assert np.abs(packet.energy - init_energy) > 1e-7


@pytest.mark.parametrize(
    "line_interaction_type",
    [
        interaction_events.LineInteractionType.SCATTER,
        interaction_events.LineInteractionType.DOWNBRANCH,
        interaction_events.LineInteractionType.MACROATOM,
    ],
)
def test_line_scatter(
    line_interaction_type,
    packet,
    verysimple_time_explosion,
    verysimple_opacity_state,
):
    init_mu = packet.mu
    init_nu = packet.nu
    init_energy = packet.energy
    full_relativity = False
    packet.initialize_line_id(
        verysimple_opacity_state, verysimple_time_explosion, full_relativity
    )
    time_explosion = verysimple_time_explosion

    line_scatter_event(
        packet,
        time_explosion,
        line_interaction_type,
        verysimple_opacity_state,
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
    verysimple_time_explosion,
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
        verysimple_time_explosion,
        full_relativity,
    )

    time_explosion = verysimple_time_explosion

    interaction_events.line_emission(
        packet,
        emission_line_id,
        time_explosion,
        verysimple_opacity_state,
        full_relativity,
    )

    assert packet.next_line_id == emission_line_id + 1
    npt.assert_almost_equal(packet.mu, expected["mu"])
    npt.assert_almost_equal(packet.energy, expected["energy"])
