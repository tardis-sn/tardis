import numpy as np
import numpy.testing as npt
import pytest
from numba import njit

from tardis.energy_input.transport.GXPacket import GXPacket, GXPacketStatus
from tardis.energy_input.transport.gamma_ray_grid import (
    calculate_distance_radial,
    distance_trace,
    move_packet,
)
from tardis.energy_input.transport.gamma_ray_interactions import scatter_type
from tardis.energy_input.util import H_CGS_KEV


@njit
def _seed_numba_random(seed: int) -> None:
    np.random.seed(seed)


def _make_gamma_packet() -> GXPacket:
    _seed_numba_random(1963)
    return GXPacket(
        location=np.array([1.0e14, 0.0, 0.0]),
        direction=np.array([1.0, 0.0, 0.0]),
        energy_rf=1.0e42,
        energy_cmf=1.0e42,
        nu_rf=1000.0 / H_CGS_KEV,
        nu_cmf=1000.0 / H_CGS_KEV,
        status=GXPacketStatus.IN_PROCESS,
        shell=0,
        time_start=1.2e5,
        time_index=0,
    )


def test_gamma_distance_radial_characterization() -> None:
    packet = _make_gamma_packet()

    distance, shell_change = calculate_distance_radial(
        packet,
        5.0e13,
        2.0e14,
    )

    npt.assert_allclose(distance, 1.0e14)
    assert shell_change == 1


def test_gamma_distance_trace_characterization() -> None:
    packet = _make_gamma_packet()

    (
        distance_interaction,
        distance_boundary,
        distance_time,
        shell_change,
    ) = distance_trace(
        packet,
        np.array([5.0e8, 1.0e9]),
        np.array([2.0e9, 3.0e9]),
        2.0e-14,
        1.0e5,
        1.5e5,
    )

    npt.assert_allclose(distance_interaction, 2.206979096801579e12)
    npt.assert_allclose(distance_boundary, 1.0e14)
    npt.assert_allclose(distance_time, 8.99377374e14)
    assert shell_change == 1


def test_gamma_move_packet_characterization() -> None:
    packet = _make_gamma_packet()

    moved_packet = move_packet(packet, 2.5e13)

    npt.assert_allclose(
        moved_packet.location,
        np.array([1.25e14, 0.0, 0.0]),
    )
    npt.assert_allclose(moved_packet.nu_cmf, 2.3339733444748378e20)
    npt.assert_allclose(moved_packet.energy_cmf, 9.652537400835258e41)


@pytest.mark.parametrize(
    ("seed", "compton_opacity", "photoabsorption_opacity", "total_opacity", "expected"),
    [
        (1963, 1.0, 0.0, 1.0, GXPacketStatus.COMPTON_SCATTER),
        (1963, 0.0, 1.0, 1.0, GXPacketStatus.PHOTOABSORPTION),
        (1963, 0.0, 0.0, 1.0, GXPacketStatus.PAIR_CREATION),
    ],
)
def test_gamma_scatter_type_characterization(
    seed: int,
    compton_opacity: float,
    photoabsorption_opacity: float,
    total_opacity: float,
    expected: GXPacketStatus,
) -> None:
    np.random.seed(seed)

    actual = scatter_type(
        compton_opacity,
        photoabsorption_opacity,
        total_opacity,
    )

    assert actual == expected
