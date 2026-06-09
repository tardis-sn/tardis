import os

import numpy as np
import numpy.testing as npt
import pytest
from numba import njit

from tardis.energy_input.transport.gamma_packet_loop import (
    gamma_packet_loop,
    process_packet_path,
)
from tardis.energy_input.transport.gamma_ray_grid import (
    calculate_distance_radial,
    distance_trace,
    move_packet,
)
from tardis.energy_input.transport.gamma_ray_interactions import (
    compton_scatter,
    get_compton_angle,
    get_compton_fraction,
    get_compton_fraction_artis,
    pair_creation_packet,
    scatter_type,
)
from tardis.energy_input.transport.GXPacket import GXPacket, GXPacketStatus
from tardis.energy_input.util import ELECTRON_MASS_ENERGY_KEV, H_CGS_KEV


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


def test_gamma_distance_radial_inward_characterization() -> None:
    packet = _make_gamma_packet()
    packet.direction = np.array([-1.0, 0.0, 0.0])

    distance, shell_change = calculate_distance_radial(
        packet,
        5.0e13,
        2.0e14,
    )

    npt.assert_allclose(distance, 5.0e13)
    # Current implementation returns +1 here because the shell-change test
    # compares against ``inner_1 or inner_2``.
    assert shell_change == 1


def test_gamma_distance_radial_no_root_characterization() -> None:
    packet = _make_gamma_packet()
    packet.location = np.array([3.0e14, 0.0, 0.0])
    packet.direction = np.array([0.0, 1.0, 0.0])

    with pytest.raises(
        ValueError, match="No root found for distance calculation!"
    ):
        calculate_distance_radial(
            packet,
            5.0e13,
            2.0e14,
        )


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


@pytest.mark.parametrize(
    (
        "tau",
        "total_opacity",
        "next_time",
        "expected_distances",
    ),
    [
        (
            1.0e-3,
            2.0e-14,
            1.5e5,
            (5.0e10, 1.0e14, 8.99377374e14),
        ),
        (
            10.0,
            2.0e-14,
            1.5e5,
            (5.0e14, 1.0e14, 8.99377374e14),
        ),
        (
            10.0,
            1.0e-20,
            1.2001e5,
            (1.0e21, 1.0e14, 2.99792458e11),
        ),
    ],
)
def test_gamma_distance_trace_winner_characterization(
    tau: float,
    total_opacity: float,
    next_time: float,
    expected_distances: tuple[float, float, float],
) -> None:
    packet = _make_gamma_packet()
    packet.tau = tau

    actual = distance_trace(
        packet,
        np.array([5.0e8, 1.0e9]),
        np.array([2.0e9, 3.0e9]),
        total_opacity,
        1.0e5,
        next_time,
    )

    npt.assert_allclose(actual[:3], expected_distances)
    assert actual[3] == 1


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


def test_gamma_compton_angle_characterization() -> None:
    _seed_numba_random(1963)

    angle, lost_energy, new_energy = get_compton_angle(1000.0)

    npt.assert_allclose(angle, 2.7743528728403617)
    npt.assert_allclose(lost_energy, 790.9444204876434)
    npt.assert_allclose(new_energy, 209.05557951235664)


def test_gamma_compton_fraction_characterization() -> None:
    _seed_numba_random(1963)

    angle, fraction = get_compton_fraction(1000.0)

    npt.assert_allclose(angle, 2.7743528728403617)
    npt.assert_allclose(fraction, 0.20905557951235665)


def test_gamma_compton_fraction_artis_characterization() -> None:
    _seed_numba_random(1963)

    angle, fraction = get_compton_fraction_artis(1000.0)

    npt.assert_allclose(angle, 2.5465507307622497)
    npt.assert_allclose(fraction, 4.577551663844115)


def test_gamma_compton_scatter_characterization() -> None:
    packet = _make_gamma_packet()
    _seed_numba_random(1963)

    actual = compton_scatter(packet, 0.5)

    npt.assert_allclose(
        actual,
        np.array([0.8710337, 0.4111476, -0.26880839]),
    )


def test_gamma_pair_creation_survives_characterization() -> None:
    packet = _make_gamma_packet()
    packet.nu_cmf = 2.0 * ELECTRON_MASS_ENERGY_KEV / H_CGS_KEV
    packet.nu_rf = packet.nu_cmf
    _seed_numba_random(2)

    actual = pair_creation_packet(packet)

    assert actual.status == GXPacketStatus.IN_PROCESS
    npt.assert_allclose(actual.nu_cmf, 1.2355899646038822e20)
    npt.assert_allclose(actual.nu_rf, 1.2261484059488753e20)
    npt.assert_allclose(actual.energy_rf, 9.923590922245192e41, rtol=1.0e-6)
    npt.assert_allclose(
        actual.direction,
        np.array([-0.27701457, -0.09836312, 0.95581778]),
    )


def test_gamma_process_packet_path_photoabsorption_characterization() -> None:
    packet = _make_gamma_packet()
    packet.status = GXPacketStatus.PHOTOABSORPTION

    actual, ejecta_energy_gained = process_packet_path(packet)

    assert actual.status == GXPacketStatus.PHOTOABSORPTION
    npt.assert_allclose(ejecta_energy_gained, 1.0e42)


def test_gamma_process_packet_path_pair_creation_characterization() -> None:
    packet = _make_gamma_packet()
    packet.status = GXPacketStatus.PAIR_CREATION
    _seed_numba_random(1963)

    actual, ejecta_energy_gained = process_packet_path(packet)

    assert actual.status == GXPacketStatus.PAIR_CREATION
    npt.assert_allclose(actual.nu_cmf, 1.2355899646038822e20)
    npt.assert_allclose(actual.energy_rf, 1.0182839191406096e42)
    npt.assert_allclose(ejecta_energy_gained, 0.0)


def test_gamma_process_packet_path_compton_characterization() -> None:
    packet = _make_gamma_packet()
    packet.status = GXPacketStatus.COMPTON_SCATTER
    _seed_numba_random(1963)

    actual, ejecta_energy_gained = process_packet_path(packet)

    assert actual.status == GXPacketStatus.PHOTOABSORPTION
    npt.assert_allclose(actual.nu_cmf, 2.4179894338175507e20)
    npt.assert_allclose(actual.energy_rf, 1.0e42)
    npt.assert_allclose(ejecta_energy_gained, 1.0e42)


def _gamma_loop_arrays() -> dict[str, np.ndarray]:
    return {
        "electron_number_density_time": np.ones((1, 2)) * 1.0e8,
        "mass_density_time": np.ones((1, 2)) * 1.0e-12,
        "iron_group_fraction_per_shell": np.array([0.5]),
        "inner_velocities": np.array([5.0e8]),
        "outer_velocities": np.array([2.0e9]),
        "dt_array": np.array([5.0e4, 5.0e4]),
        "times": np.array([1.0e5, 1.5e5, 2.0e5]),
        "effective_time_array": np.array([1.25e5, 1.75e5]),
        "energy_bins": np.array([0.0, 500.0, 1500.0, 3000.0]),
        "energy_out": np.zeros((3, 2)),
        "energy_out_cosi": np.zeros((3, 2)),
        "total_energy": np.zeros((1, 2)),
        "energy_deposited_gamma": np.zeros((1, 2)),
        "packets_info_array": np.zeros((1, 8)),
    }


def test_gamma_packet_loop_negative_time_index_characterization() -> None:
    packet = _make_gamma_packet()
    packet.time_index = -1

    with pytest.raises(ValueError, match="Packet time index less than 0!"):
        gamma_packet_loop(
            [packet],
            -1.0,
            "tardis",
            "artis",
            **_gamma_loop_arrays(),
        )


@pytest.mark.parametrize("grey_opacity", [-1.0, 0.1])
def test_gamma_packet_loop_escape_binning_characterization(
    grey_opacity: float,
) -> None:
    if os.environ.get("NUMBA_DISABLE_JIT") == "1" and grey_opacity >= 0.0:
        pytest.xfail(
            "Current Python execution path leaves doppler_factor undefined for grey opacity."
        )

    packet = _make_gamma_packet()
    packet.location = np.array([1.9e14, 0.0, 0.0])
    packet.direction = np.array([1.0, 0.0, 0.0])
    arrays = _gamma_loop_arrays()
    _seed_numba_random(1963)

    (
        energy_out,
        energy_out_cosi,
        packets_info_array,
        energy_deposited_gamma,
        total_energy,
    ) = gamma_packet_loop(
        [packet],
        grey_opacity,
        "kasen",
        "artis",
        **arrays,
    )

    assert packet.status == GXPacketStatus.ESCAPED
    assert packet.shell == 1
    npt.assert_allclose(energy_out[1, 0], 8.27133474e16)
    npt.assert_allclose(energy_out_cosi[1, 0], 2.0e-8)
    npt.assert_allclose(packets_info_array[0, 1], GXPacketStatus.ESCAPED)
    npt.assert_allclose(packets_info_array[0, 5], 2.0e37)
    npt.assert_allclose(packets_info_array[0, 6], 1.0e42)
    npt.assert_allclose(energy_deposited_gamma, np.zeros((1, 2)))
    npt.assert_allclose(total_energy, np.zeros((1, 2)))
