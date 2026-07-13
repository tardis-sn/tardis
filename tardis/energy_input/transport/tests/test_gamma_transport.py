import os

import numpy as np
import pytest

import tardis.energy_input.transport.gamma_packet_loop as gamma_loop_module
from tardis.conftest import sync_ndarray_assert_allclose
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

RTOL = 1.0e-12


@pytest.fixture
def gamma_packet(basic_gamma_ray: GXPacket, set_seed_fixture) -> GXPacket:
    set_seed_fixture(1963)
    basic_gamma_ray.location = np.array([1.0e14, 0.0, 0.0])
    basic_gamma_ray.direction = np.array([1.0, 0.0, 0.0])
    basic_gamma_ray.energy_rf = 1.0e42
    basic_gamma_ray.energy_cmf = 1.0e42
    basic_gamma_ray.nu_rf = 1000.0 / H_CGS_KEV
    basic_gamma_ray.nu_cmf = 1000.0 / H_CGS_KEV
    basic_gamma_ray.status = GXPacketStatus.IN_PROCESS
    basic_gamma_ray.shell = 0
    basic_gamma_ray.time_start = 1.2e5
    basic_gamma_ray.time_idx = 0
    return basic_gamma_ray


@pytest.fixture
def gamma_loop_arrays() -> dict[str, np.ndarray]:
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


def test_gamma_distance_radial(
    gamma_packet: GXPacket,
    regression_data,
) -> None:
    packet = gamma_packet

    distance, shell_change = calculate_distance_radial(
        packet,
        5.0e13,
        2.0e14,
    )

    sync_ndarray_assert_allclose(regression_data, distance, rtol=RTOL)
    assert shell_change == 1


def test_gamma_distance_radial_inward(
    gamma_packet: GXPacket,
    regression_data,
) -> None:
    packet = gamma_packet
    packet.direction = np.array([-1.0, 0.0, 0.0])

    distance, shell_change = calculate_distance_radial(
        packet,
        5.0e13,
        2.0e14,
    )

    sync_ndarray_assert_allclose(regression_data, distance, rtol=RTOL)
    # Current implementation returns +1 here because the shell-change test
    # compares against ``inner_1 or inner_2``.
    assert shell_change == 1


def test_gamma_distance_radial_no_root(
    gamma_packet: GXPacket,
) -> None:
    packet = gamma_packet
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


def test_gamma_distance_trace(
    gamma_packet: GXPacket,
    regression_data,
) -> None:
    (
        distance_interaction,
        distance_boundary,
        distance_time,
        shell_change,
    ) = distance_trace(
        gamma_packet,
        np.array([5.0e8, 1.0e9]),
        np.array([2.0e9, 3.0e9]),
        2.0e-14,
        1.0e5,
        1.5e5,
    )

    sync_ndarray_assert_allclose(
        regression_data,
        np.array([distance_interaction, distance_boundary, distance_time]),
        rtol=RTOL,
    )
    assert shell_change == 1


@pytest.mark.parametrize(
    (
        "tau",
        "total_opacity",
        "next_time",
    ),
    [
        (1.0e-3, 2.0e-14, 1.5e5),
        (10.0, 2.0e-14, 1.5e5),
        (10.0, 1.0e-20, 1.2001e5),
    ],
)
def test_gamma_distance_trace_winner(
    gamma_packet: GXPacket,
    tau: float,
    total_opacity: float,
    next_time: float,
    regression_data,
) -> None:
    gamma_packet.tau = tau

    actual = distance_trace(
        gamma_packet,
        np.array([5.0e8, 1.0e9]),
        np.array([2.0e9, 3.0e9]),
        total_opacity,
        1.0e5,
        next_time,
    )

    sync_ndarray_assert_allclose(regression_data, actual[:3], rtol=RTOL)
    assert actual[3] == 1


def test_gamma_move_packet(gamma_packet: GXPacket, regression_data) -> None:
    moved_packet = move_packet(gamma_packet, 2.5e13)

    sync_ndarray_assert_allclose(
        regression_data,
        moved_packet.location,
        moved_packet.nu_cmf,
        moved_packet.energy_cmf,
        rtol=RTOL,
    )


@pytest.mark.parametrize(
    (
        "seed",
        "compton_opacity",
        "photoabsorption_opacity",
        "total_opacity",
        "expected",
    ),
    [
        (1963, 1.0, 0.0, 1.0, GXPacketStatus.COMPTON_SCATTER),
        (1963, 0.0, 1.0, 1.0, GXPacketStatus.PHOTOABSORPTION),
        (1963, 0.0, 0.0, 1.0, GXPacketStatus.PAIR_CREATION),
    ],
)
def test_gamma_scatter_type(
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


def test_gamma_compton_angle(set_seed_fixture, regression_data) -> None:
    set_seed_fixture(1963)

    angle, lost_energy, new_energy = get_compton_angle(1000.0)

    sync_ndarray_assert_allclose(
        regression_data, np.array([angle, lost_energy, new_energy]), rtol=RTOL
    )


def test_gamma_compton_fraction(set_seed_fixture, regression_data) -> None:
    set_seed_fixture(1963)

    angle, fraction = get_compton_fraction(1000.0)

    sync_ndarray_assert_allclose(
        regression_data, np.array([angle, fraction]), rtol=RTOL
    )


def test_gamma_compton_fraction_artis(
    set_seed_fixture, regression_data
) -> None:
    set_seed_fixture(1963)

    angle, fraction = get_compton_fraction_artis(1000.0)

    sync_ndarray_assert_allclose(
        regression_data, np.array([angle, fraction]), rtol=RTOL
    )


def test_gamma_compton_scatter(
    gamma_packet: GXPacket,
    set_seed_fixture,
    regression_data,
) -> None:
    set_seed_fixture(1963)

    actual = compton_scatter(gamma_packet, 0.5)

    sync_ndarray_assert_allclose(regression_data, actual, rtol=RTOL)


def test_gamma_pair_creation_survives(
    gamma_packet: GXPacket,
    set_seed_fixture,
    regression_data,
) -> None:
    packet = gamma_packet
    packet.nu_cmf = 2.0 * ELECTRON_MASS_ENERGY_KEV / H_CGS_KEV
    packet.nu_rf = packet.nu_cmf
    set_seed_fixture(2)

    actual = pair_creation_packet(packet)

    assert actual.status == GXPacketStatus.IN_PROCESS
    sync_ndarray_assert_allclose(
        regression_data,
        actual.nu_cmf,
        actual.nu_rf,
        actual.energy_rf,
        actual.direction,
        rtol=RTOL,
    )


def test_gamma_process_packet_path_photoabsorption(
    gamma_packet: GXPacket,
    regression_data,
) -> None:
    packet = gamma_packet
    packet.status = GXPacketStatus.PHOTOABSORPTION

    actual, ejecta_energy_gained = process_packet_path(packet)

    assert actual.status == GXPacketStatus.PHOTOABSORPTION
    sync_ndarray_assert_allclose(
        regression_data, ejecta_energy_gained, rtol=RTOL
    )


def test_gamma_process_packet_path_pair_creation(
    gamma_packet: GXPacket,
    set_seed_fixture,
    regression_data,
) -> None:
    packet = gamma_packet
    packet.status = GXPacketStatus.PAIR_CREATION
    set_seed_fixture(1963)

    actual, ejecta_energy_gained = process_packet_path(packet)

    assert actual.status == GXPacketStatus.PAIR_CREATION
    sync_ndarray_assert_allclose(
        regression_data,
        actual.nu_cmf,
        actual.energy_rf,
        ejecta_energy_gained,
        rtol=RTOL,
    )


def test_gamma_process_packet_path_compton(
    gamma_packet: GXPacket,
    set_seed_fixture,
    regression_data,
) -> None:
    packet = gamma_packet
    packet.status = GXPacketStatus.COMPTON_SCATTER
    set_seed_fixture(1963)

    actual, ejecta_energy_gained = process_packet_path(packet)

    assert actual.status == GXPacketStatus.PHOTOABSORPTION
    sync_ndarray_assert_allclose(
        regression_data,
        actual.nu_cmf,
        actual.energy_rf,
        ejecta_energy_gained,
        rtol=RTOL,
    )


def test_gamma_packet_loop_negative_time_index(
    gamma_packet: GXPacket,
    gamma_loop_arrays: dict[str, np.ndarray],
) -> None:
    gamma_packet.time_idx = -1

    with pytest.raises(ValueError, match="Packet time index less than 0!"):
        gamma_packet_loop(
            [gamma_packet],
            -1.0,
            "tardis",
            "artis",
            **gamma_loop_arrays,
        )


@pytest.mark.parametrize("grey_opacity", [-1.0, 0.1])
def test_gamma_packet_loop_escape_binning(
    gamma_packet: GXPacket,
    gamma_loop_arrays: dict[str, np.ndarray],
    grey_opacity: float,
    set_seed_fixture,
    regression_data,
) -> None:
    # Put the packet just inside the outer boundary so the loop exercises
    # escape binning rather than an interaction branch.
    if os.environ.get("NUMBA_DISABLE_JIT") == "1" and grey_opacity >= 0.0:
        pytest.xfail(
            "Current Python execution path leaves doppler_factor undefined "
            "for grey opacity."
        )

    packet = gamma_packet
    packet.location = np.array([1.9e14, 0.0, 0.0])
    packet.direction = np.array([1.0, 0.0, 0.0])
    set_seed_fixture(1963)

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
        **gamma_loop_arrays,
    )

    assert packet.status == GXPacketStatus.ESCAPED
    assert packet.shell == 1
    assert packets_info_array[0, 1] == GXPacketStatus.ESCAPED
    sync_ndarray_assert_allclose(
        regression_data,
        energy_out[1, 0],
        energy_out_cosi[1, 0],
        packets_info_array[0, 5],
        packets_info_array[0, 6],
        energy_deposited_gamma,
        total_energy,
        rtol=RTOL,
    )


def test_gamma_packet_loop_tardis_opacity(
    gamma_packet: GXPacket,
    gamma_loop_arrays: dict[str, np.ndarray],
    set_seed_fixture,
    regression_data,
) -> None:
    # Put the packet just inside the outer boundary so the TARDIS opacity path
    # reaches escape binning deterministically.
    packet = gamma_packet
    packet.location = np.array([1.9e14, 0.0, 0.0])
    packet.direction = np.array([1.0, 0.0, 0.0])
    set_seed_fixture(1963)

    (
        energy_out,
        energy_out_cosi,
        packets_info_array,
        energy_deposited_gamma,
        total_energy,
    ) = gamma_packet_loop(
        [packet],
        -1.0,
        "tardis",
        "tardis",
        **gamma_loop_arrays,
    )

    assert packet.status == GXPacketStatus.ESCAPED
    assert packets_info_array[0, 1] == GXPacketStatus.ESCAPED
    sync_ndarray_assert_allclose(
        regression_data,
        energy_out[1, 0],
        energy_out_cosi[1, 0],
        energy_deposited_gamma,
        total_energy,
        rtol=RTOL,
    )


@pytest.mark.parametrize(
    ("photoabsorption_opacity_type", "pair_creation_opacity_type", "match"),
    [
        ("invalid", "artis", "Invalid photoabsorption opacity type!"),
        ("kasen", "invalid", "Invalid pair creation opacity type!"),
    ],
)
def test_gamma_packet_loop_invalid_opacity_type(
    gamma_packet: GXPacket,
    gamma_loop_arrays: dict[str, np.ndarray],
    photoabsorption_opacity_type: str,
    pair_creation_opacity_type: str,
    match: str,
) -> None:
    with pytest.raises(ValueError, match=match):
        gamma_packet_loop(
            [gamma_packet],
            -1.0,
            photoabsorption_opacity_type,
            pair_creation_opacity_type,
            **gamma_loop_arrays,
        )


def test_gamma_packet_loop_time_boundary_end_numba_disabled(
    monkeypatch,
    python_numba_disabled,
    gamma_packet: GXPacket,
    gamma_loop_arrays: dict[str, np.ndarray],
    regression_data,
) -> None:
    packet = gamma_packet
    packet.time_idx = 1

    # Force the time-boundary branch without depending on sampled distances.
    monkeypatch.setattr(
        gamma_loop_module,
        "distance_trace",
        lambda *args: (10.0, 20.0, 1.0, 0),
    )

    with pytest.raises(UnboundLocalError):
        gamma_packet_loop(
            [packet],
            -1.0,
            "kasen",
            "artis",
            **gamma_loop_arrays,
        )

    assert packet.status == GXPacketStatus.END
    assert packet.shell == 0
    sync_ndarray_assert_allclose(
        regression_data,
        gamma_loop_arrays["energy_deposited_gamma"],
        gamma_loop_arrays["total_energy"],
        rtol=RTOL,
    )


def test_gamma_packet_loop_inner_boundary_end_numba_disabled(
    monkeypatch,
    python_numba_disabled,
    gamma_packet: GXPacket,
    gamma_loop_arrays: dict[str, np.ndarray],
    regression_data,
) -> None:
    packet = gamma_packet

    # Force the inner-boundary branch without depending on sampled distances.
    monkeypatch.setattr(
        gamma_loop_module,
        "distance_trace",
        lambda *args: (20.0, 1.0, 10.0, -1),
    )

    with pytest.raises(UnboundLocalError):
        gamma_packet_loop(
            [packet],
            -1.0,
            "kasen",
            "artis",
            **gamma_loop_arrays,
        )

    assert packet.status == GXPacketStatus.END
    assert packet.shell == -1
    sync_ndarray_assert_allclose(
        regression_data,
        packet.energy_rf,
        packet.energy_cmf,
        gamma_loop_arrays["energy_deposited_gamma"],
        gamma_loop_arrays["total_energy"],
        rtol=RTOL,
    )


def test_gamma_packet_loop_interaction_deposition_numba_disabled(
    monkeypatch,
    python_numba_disabled,
    gamma_packet: GXPacket,
    gamma_loop_arrays: dict[str, np.ndarray],
    regression_data,
) -> None:
    packet = gamma_packet

    # Force an immediate photoabsorption interaction to characterize deposition
    # bookkeeping in the loop.
    monkeypatch.setattr(
        gamma_loop_module,
        "distance_trace",
        lambda *args: (1.0, 20.0, 10.0, 0),
    )
    monkeypatch.setattr(
        gamma_loop_module,
        "scatter_type",
        lambda *args: GXPacketStatus.PHOTOABSORPTION,
    )

    _, _, packets_info_array, energy_deposited_gamma, total_energy = (
        gamma_packet_loop(
            [packet],
            -1.0,
            "kasen",
            "artis",
            **gamma_loop_arrays,
        )
    )

    assert packet.status == GXPacketStatus.PHOTOABSORPTION
    sync_ndarray_assert_allclose(
        regression_data,
        energy_deposited_gamma,
        total_energy,
        packets_info_array,
        rtol=RTOL,
    )


def test_gamma_packet_loop_scattered_escape_numba_disabled(
    monkeypatch,
    python_numba_disabled,
    gamma_packet: GXPacket,
    gamma_loop_arrays: dict[str, np.ndarray],
    regression_data,
) -> None:
    packet = gamma_packet
    distances_to_return = [
        # First loop step: interaction distance wins.
        (1.0, 20.0, 10.0, 0),
        # Second loop step: boundary distance wins and moves out one shell.
        (20.0, 1.0, 10.0, 1),
    ]

    def distance_trace_interacts_then_escapes(
        *args: object,
    ) -> tuple[float, float, float, int]:
        assert distances_to_return
        return distances_to_return.pop(0)

    def process_packet_path_as_pair_creation(
        packet: GXPacket,
    ) -> tuple[GXPacket, float]:
        packet.status = GXPacketStatus.PAIR_CREATION
        return packet, 0.0

    # Force an interaction followed by boundary escape to characterize the
    # status reset after a scattering-like path.
    monkeypatch.setattr(
        gamma_loop_module,
        "distance_trace",
        distance_trace_interacts_then_escapes,
    )
    monkeypatch.setattr(
        gamma_loop_module,
        "scatter_type",
        lambda *args: GXPacketStatus.PAIR_CREATION,
    )
    monkeypatch.setattr(
        gamma_loop_module,
        "process_packet_path",
        process_packet_path_as_pair_creation,
    )

    with pytest.raises(UnboundLocalError):
        gamma_packet_loop(
            [packet],
            -1.0,
            "kasen",
            "artis",
            **gamma_loop_arrays,
        )

    assert packet.status == GXPacketStatus.IN_PROCESS
    sync_ndarray_assert_allclose(
        regression_data,
        gamma_loop_arrays["energy_deposited_gamma"],
        gamma_loop_arrays["total_energy"],
        rtol=RTOL,
    )
