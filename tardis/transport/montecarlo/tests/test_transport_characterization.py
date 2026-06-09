import numpy as np
import numpy.testing as npt
import pytest
from numba import njit

from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.model.geometry.radial1d_nonhomologous import (
    NumbaNonhomologousRadial1DGeometry,
)
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.opacities.opacity_state_numba_iip import OpacityStateNumbaIIP
from tardis.transport.montecarlo import RPacket
from tardis.transport.montecarlo.estimators.estimators_bulk import (
    init_estimators_bulk,
)
from tardis.transport.montecarlo.estimators.estimators_line import (
    EstimatorsLine,
)
from tardis.transport.montecarlo.modes.classic.rad_packet_transport import (
    move_packet_across_shell_boundary as classic_move_boundary,
)
from tardis.transport.montecarlo.modes.classic.rad_packet_transport import (
    move_r_packet as classic_move_r_packet,
)
from tardis.transport.montecarlo.modes.classic.rad_packet_transport import (
    trace_packet as classic_trace_packet,
)
from tardis.transport.montecarlo.modes.iip.rad_packet_transport import (
    move_packet_across_shell_boundary as iip_move_boundary,
)
from tardis.transport.montecarlo.modes.iip.rad_packet_transport import (
    move_r_packet as iip_move_r_packet,
)
from tardis.transport.montecarlo.modes.iip.rad_packet_transport import (
    trace_packet as iip_trace_packet,
)
from tardis.transport.montecarlo.modes.nonhomologous.rad_packet_transport import (
    move_packet_across_shell_boundary as nonhomologous_move_boundary,
)
from tardis.transport.montecarlo.modes.nonhomologous.rad_packet_transport import (
    move_r_packet as nonhomologous_move_r_packet,
)
from tardis.transport.montecarlo.modes.nonhomologous.rad_packet_transport import (
    trace_packet as nonhomologous_trace_packet,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
    PacketStatus,
)


@njit
def _seed_numba_random(seed: int) -> None:
    np.random.seed(seed)


@pytest.fixture
def r_packet() -> RPacket:
    packet = RPacket(
        r=7.5e14,
        mu=0.3,
        nu=4.0e14,
        energy=0.9,
        seed=1963,
        index=0,
    )
    packet.current_shell_id = 1
    return packet


def _make_rpacket(
    *,
    current_shell_id: int = 0,
    next_line_id: int = 0,
    prev_line_id: int = 0,
) -> RPacket:
    packet = RPacket(
        r=7.5e14,
        mu=0.3,
        nu=4.0e14,
        energy=0.9,
        seed=1963,
        index=0,
    )
    packet.current_shell_id = current_shell_id
    packet.next_line_id = next_line_id
    packet.prev_line_id = prev_line_id
    return packet


def _make_radial_geometry(r_outer_first_shell: float) -> NumbaRadial1DGeometry:
    return NumbaRadial1DGeometry(
        np.array([7.0e14, 8.0e14]),
        np.array([r_outer_first_shell, 3.0e16]),
        np.array([-1.0, -1.0]),
        np.array([-1.0, -1.0]),
    )


def _make_nonhomologous_geometry(
    *,
    negative_velocity_gradient: bool = False,
    r_outer_first_shell: float = 8.0e14,
) -> NumbaNonhomologousRadial1DGeometry:
    if negative_velocity_gradient:
        v_inner = np.array([1.5e9, 2.0e9])
        v_outer = np.array([1.0e9, 1.5e9])
    else:
        v_inner = np.array([1.0e9, 1.5e9])
        v_outer = np.array([1.5e9, 2.0e9])

    return NumbaNonhomologousRadial1DGeometry(
        np.array([7.0e14, 8.0e14]),
        np.array([r_outer_first_shell, 3.0e16]),
        v_inner,
        v_outer,
    )


def _make_line_estimators(
    no_of_lines: int = 2, no_of_shells: int = 2
) -> EstimatorsLine:
    return EstimatorsLine(
        np.zeros((no_of_lines, no_of_shells)),
        np.zeros((no_of_lines, no_of_shells)),
    )


def _opacity_state_args(
    line_list_nu: np.ndarray, tau_sobolev: np.ndarray
) -> tuple:
    no_of_lines = len(line_list_nu)
    no_of_shells = tau_sobolev.shape[1]
    return (
        np.ones(no_of_shells) * 1.0e8,
        np.ones(no_of_shells) * 1.0e4,
        line_list_nu,
        tau_sobolev,
        np.zeros((1, no_of_shells)),
        np.zeros(no_of_lines, dtype=np.int64),
        np.zeros(1, dtype=np.int64),
        np.zeros(1, dtype=np.int64),
        np.zeros(1, dtype=np.int64),
        np.zeros(1, dtype=np.int64),
        np.zeros(1),
        np.zeros((1, no_of_shells)),
        np.zeros(1),
        np.zeros(1),
        np.zeros(1, dtype=np.int64),
        np.zeros((1, no_of_shells)),
        np.zeros(1),
        np.zeros(1),
        np.zeros(no_of_shells),
        np.zeros((1, no_of_shells)),
        np.zeros(1, dtype=np.int64),
        0,
    )


def _make_opacity_state(
    line_list_nu: list[float],
    tau_sobolev: np.ndarray,
) -> OpacityStateNumba:
    return OpacityStateNumba(
        *_opacity_state_args(np.array(line_list_nu), tau_sobolev)
    )


def _make_iip_opacity_state(
    line_list_nu: list[float],
    tau_sobolev: np.ndarray,
) -> OpacityStateNumbaIIP:
    no_of_shells = tau_sobolev.shape[1]
    return OpacityStateNumbaIIP(
        *_opacity_state_args(np.array(line_list_nu), tau_sobolev),
        np.ones((no_of_shells, 1, 1)),
    )


@pytest.mark.parametrize(
    "move_r_packet",
    [classic_move_r_packet, iip_move_r_packet],
)
@pytest.mark.parametrize(
    (
        "enable_full_relativity",
        "expected_mean_intensity_total",
        "expected_mean_frequency",
    ),
    [
        (
            False,
            np.array([0.0, 8.998701024436969e12, 0.0]),
            np.array([0.0, 3.5989608945423537e27, 0.0]),
        ),
        (
            True,
            np.array([0.0, 8.997404318887818e12, 0.0]),
            np.array([0.0, 3.598442703630763e27, 0.0]),
        ),
    ],
)
def test_homologous_move_r_packet_characterization(
    move_r_packet,
    r_packet: RPacket,
    enable_full_relativity: bool,
    expected_mean_intensity_total: np.ndarray,
    expected_mean_frequency: np.ndarray,
) -> None:
    packet = r_packet
    estimators = init_estimators_bulk(3)

    move_r_packet(
        packet,
        1.0e13,
        5.2e7,
        estimators,
        enable_full_relativity,
    )

    npt.assert_allclose(packet.r, 7.530604225425739e14)
    npt.assert_allclose(packet.mu, 0.3120599529139568)
    npt.assert_allclose(
        estimators.mean_intensity_total,
        expected_mean_intensity_total,
    )
    npt.assert_allclose(
        estimators.mean_frequency,
        expected_mean_frequency,
    )


def test_nonhomologous_move_r_packet_characterization(
    r_packet: RPacket,
) -> None:
    geometry = NumbaNonhomologousRadial1DGeometry(
        np.array([7.0e14, 8.0e14]),
        np.array([8.0e14, 9.0e14]),
        np.array([1.0e9, 1.5e9]),
        np.array([1.5e9, 2.0e9]),
    )
    packet = r_packet
    packet.current_shell_id = 0
    estimators = init_estimators_bulk(2)

    nonhomologous_move_r_packet(
        packet,
        1.0e13,
        geometry,
        estimators,
        False,
    )

    npt.assert_allclose(packet.r, 7.530604225425739e14)
    npt.assert_allclose(packet.mu, 0.3120599529139568)
    npt.assert_allclose(
        estimators.mean_intensity_total,
        np.array([8.887422117870625e12, 0.0]),
    )
    npt.assert_allclose(
        estimators.mean_frequency,
        np.array([3.510500973387377e27, 0.0]),
    )


def test_nonhomologous_move_r_packet_full_relativity_characterization(
    r_packet: RPacket,
) -> None:
    geometry = NumbaNonhomologousRadial1DGeometry(
        np.array([7.0e14, 8.0e14]),
        np.array([8.0e14, 9.0e14]),
        np.array([1.0e9, 1.5e9]),
        np.array([1.5e9, 2.0e9]),
    )
    packet = r_packet
    packet.current_shell_id = 0
    estimators = init_estimators_bulk(2)

    with pytest.raises(
        NotImplementedError,
        match=r"Full relativity not implemented for non-homologous mode.",
    ):
        nonhomologous_move_r_packet(
            packet,
            1.0e13,
            geometry,
            estimators,
            True,
        )


@pytest.mark.parametrize(
    "move_r_packet",
    [classic_move_r_packet, iip_move_r_packet],
)
def test_homologous_move_r_packet_zero_distance_characterization(
    move_r_packet,
    r_packet: RPacket,
) -> None:
    packet = r_packet
    estimators = init_estimators_bulk(3)

    move_r_packet(
        packet,
        0.0,
        5.2e7,
        estimators,
        False,
    )

    npt.assert_allclose(packet.r, 7.5e14)
    npt.assert_allclose(packet.mu, 0.3)
    npt.assert_allclose(estimators.mean_intensity_total, np.zeros(3))
    npt.assert_allclose(estimators.mean_frequency, np.zeros(3))


def test_nonhomologous_move_r_packet_zero_distance_characterization(
    r_packet: RPacket,
) -> None:
    geometry = _make_nonhomologous_geometry()
    packet = r_packet
    packet.current_shell_id = 0
    estimators = init_estimators_bulk(2)

    nonhomologous_move_r_packet(
        packet,
        0.0,
        geometry,
        estimators,
        False,
    )

    npt.assert_allclose(packet.r, 7.5e14)
    npt.assert_allclose(packet.mu, 0.3)
    npt.assert_allclose(estimators.mean_intensity_total, np.zeros(2))
    npt.assert_allclose(estimators.mean_frequency, np.zeros(2))


@pytest.mark.parametrize(
    "move_boundary",
    [classic_move_boundary, iip_move_boundary, nonhomologous_move_boundary],
)
@pytest.mark.parametrize(
    (
        "current_shell_id",
        "delta_shell",
        "no_of_shells",
        "expected_status",
        "expected_shell_id",
    ),
    [
        (2, 1, 3, PacketStatus.EMITTED, 2),
        (0, -1, 3, PacketStatus.REABSORBED, 0),
        (1, 0, 3, PacketStatus.IN_PROCESS, 1),
        (1, 1, 3, PacketStatus.IN_PROCESS, 2),
    ],
)
def test_move_packet_across_shell_boundary_characterization(
    move_boundary,
    r_packet: RPacket,
    current_shell_id: int,
    delta_shell: int,
    no_of_shells: int,
    expected_status: PacketStatus,
    expected_shell_id: int,
) -> None:
    packet = r_packet
    packet.current_shell_id = current_shell_id

    move_boundary(packet, delta_shell, no_of_shells)

    assert packet.status == expected_status
    assert packet.current_shell_id == expected_shell_id


@pytest.mark.parametrize(
    (
        "opacity_electron",
        "tau_sobolev",
        "disable_line_scattering",
        "r_outer_first_shell",
        "line_list_nu",
        "expected_distance",
        "expected_interaction_type",
        "expected_j_blue",
        "expected_edot",
    ),
    [
        (
            1.0e-20,
            np.zeros((2, 2)),
            False,
            8.0e14,
            [3.95e14, 3.90e14],
            1.3294552658190883e14,
            InteractionType.BOUNDARY,
            np.zeros((2, 2)),
            np.zeros((2, 2)),
        ),
        (
            1.0e-12,
            np.zeros((2, 2)),
            False,
            8.0e14,
            [3.95e14, 3.90e14],
            4.4139581936031586e10,
            InteractionType.ESCATTERING,
            np.zeros((2, 2)),
            np.zeros((2, 2)),
        ),
        (
            1.0e-20,
            np.ones((2, 2)) * 100.0,
            False,
            2.0e16,
            [3.999e14, 3.998e14],
            1.647301954000838e14,
            InteractionType.LINE,
            np.array([[2.2494375e-15, 0.0], [0.0, 0.0]]),
            np.array([[0.899775, 0.0], [0.0, 0.0]]),
        ),
        (
            1.0e-20,
            np.ones((2, 2)) * 100.0,
            True,
            2.0e16,
            [3.999e14, 3.998e14],
            -9.995586041806397e21,
            InteractionType.ESCATTERING,
            np.array([[2.2494375e-15, 0.0], [0.0, 0.0]]),
            np.array([[0.899775, 0.0], [0.0, 0.0]]),
        ),
    ],
)
def test_classic_trace_packet_characterization(
    opacity_electron: float,
    tau_sobolev: np.ndarray,
    disable_line_scattering: bool,
    r_outer_first_shell: float,
    line_list_nu: list[float],
    expected_distance: float,
    expected_interaction_type: InteractionType,
    expected_j_blue: np.ndarray,
    expected_edot: np.ndarray,
) -> None:
    packet = _make_rpacket()
    geometry = _make_radial_geometry(r_outer_first_shell)
    opacity_state = _make_opacity_state(line_list_nu, tau_sobolev)
    estimators = _make_line_estimators()
    _seed_numba_random(1963)

    distance, interaction_type, delta_shell = classic_trace_packet(
        packet,
        geometry,
        5.2e7,
        opacity_state,
        estimators,
        opacity_electron,
        False,
        disable_line_scattering,
    )

    npt.assert_allclose(distance, expected_distance)
    assert interaction_type == expected_interaction_type
    assert delta_shell == 1
    assert packet.next_line_id == (1 if disable_line_scattering else 0)
    npt.assert_allclose(estimators.mean_intensity_blueward, expected_j_blue)
    npt.assert_allclose(estimators.energy_deposition_line_rate, expected_edot)


def test_classic_trace_packet_no_line_fallthrough_characterization() -> None:
    packet = _make_rpacket(next_line_id=2)
    geometry = _make_radial_geometry(2.0e16)
    opacity_state = _make_opacity_state(
        [3.999e14, 3.998e14], np.zeros((2, 2))
    )
    estimators = _make_line_estimators()
    _seed_numba_random(1963)

    distance, interaction_type, delta_shell = classic_trace_packet(
        packet,
        geometry,
        5.2e7,
        opacity_state,
        estimators,
        1.0e-12,
        False,
        False,
    )

    npt.assert_allclose(distance, 4.4139581936031586e10)
    assert interaction_type == InteractionType.ESCATTERING
    assert delta_shell == 1
    assert packet.next_line_id == 2
    npt.assert_allclose(estimators.mean_intensity_blueward, np.zeros((2, 2)))


@pytest.mark.parametrize(
    (
        "chi_continuum",
        "escat_prob",
        "tau_sobolev",
        "disable_line_scattering",
        "r_outer_first_shell",
        "line_list_nu",
        "expected_distance",
        "expected_interaction_type",
        "expected_j_blue",
        "expected_edot",
    ),
    [
        (
            1.0e-20,
            0.5,
            np.zeros((2, 2)),
            False,
            8.0e14,
            [3.95e14, 3.90e14],
            1.3294552658190883e14,
            InteractionType.BOUNDARY,
            np.zeros((2, 2)),
            np.zeros((2, 2)),
        ),
        (
            1.0e-12,
            1.0,
            np.zeros((2, 2)),
            False,
            8.0e14,
            [3.95e14, 3.90e14],
            4.4139581936031586e10,
            InteractionType.ESCATTERING,
            np.zeros((2, 2)),
            np.zeros((2, 2)),
        ),
        (
            1.0e-12,
            0.0,
            np.zeros((2, 2)),
            False,
            8.0e14,
            [3.95e14, 3.90e14],
            4.4139581936031586e10,
            InteractionType.CONTINUUM_PROCESS,
            np.zeros((2, 2)),
            np.zeros((2, 2)),
        ),
        (
            1.0e-20,
            0.5,
            np.ones((2, 2)) * 100.0,
            False,
            2.0e16,
            [3.999e14, 3.998e14],
            1.647301954000838e14,
            InteractionType.LINE,
            np.array([[2.2494375e-15, 0.0], [0.0, 0.0]]),
            np.array([[0.899775, 0.0], [0.0, 0.0]]),
        ),
        (
            1.0e-20,
            0.5,
            np.ones((2, 2)) * 100.0,
            True,
            2.0e16,
            [3.999e14, 3.998e14],
            -9.995586041806397e21,
            InteractionType.ESCATTERING,
            np.array([[2.2494375e-15, 0.0], [0.0, 0.0]]),
            np.array([[0.899775, 0.0], [0.0, 0.0]]),
        ),
    ],
)
def test_iip_trace_packet_characterization(
    chi_continuum: float,
    escat_prob: float,
    tau_sobolev: np.ndarray,
    disable_line_scattering: bool,
    r_outer_first_shell: float,
    line_list_nu: list[float],
    expected_distance: float,
    expected_interaction_type: InteractionType,
    expected_j_blue: np.ndarray,
    expected_edot: np.ndarray,
) -> None:
    packet = _make_rpacket()
    geometry = _make_radial_geometry(r_outer_first_shell)
    opacity_state = _make_iip_opacity_state(line_list_nu, tau_sobolev)
    estimators = _make_line_estimators()
    _seed_numba_random(1963)

    distance, interaction_type, delta_shell = iip_trace_packet(
        packet,
        geometry,
        5.2e7,
        opacity_state,
        estimators,
        chi_continuum,
        escat_prob,
        False,
        disable_line_scattering,
    )

    npt.assert_allclose(distance, expected_distance)
    assert interaction_type == expected_interaction_type
    assert delta_shell == 1
    assert packet.next_line_id == (1 if disable_line_scattering else 0)
    npt.assert_allclose(estimators.mean_intensity_blueward, expected_j_blue)
    npt.assert_allclose(estimators.energy_deposition_line_rate, expected_edot)


def test_iip_trace_packet_no_line_fallthrough_characterization() -> None:
    packet = _make_rpacket(next_line_id=2)
    geometry = _make_radial_geometry(2.0e16)
    opacity_state = _make_iip_opacity_state(
        [3.999e14, 3.998e14], np.zeros((2, 2))
    )
    estimators = _make_line_estimators()
    _seed_numba_random(1963)

    distance, interaction_type, delta_shell = iip_trace_packet(
        packet,
        geometry,
        5.2e7,
        opacity_state,
        estimators,
        1.0e-12,
        1.0,
        False,
        False,
    )

    npt.assert_allclose(distance, 4.4139581936031586e10)
    assert interaction_type == InteractionType.ESCATTERING
    assert delta_shell == 1
    assert packet.next_line_id == 2
    npt.assert_allclose(estimators.mean_intensity_blueward, np.zeros((2, 2)))


@pytest.mark.parametrize(
    (
        "geometry",
        "opacity_electron",
        "tau_sobolev",
        "line_list_nu",
        "next_line_id",
        "prev_line_id",
        "expected_distance",
        "expected_interaction_type",
        "expected_next_line_id",
        "expected_prev_line_id",
        "expected_j_blue",
        "expected_edot",
    ),
    [
        (
            _make_nonhomologous_geometry(),
            1.0e-20,
            np.zeros((2, 2)),
            [3.95e14, 3.90e14],
            0,
            0,
            1.3294552658190883e14,
            InteractionType.BOUNDARY,
            0,
            -1,
            np.zeros((2, 2)),
            np.zeros((2, 2)),
        ),
        (
            _make_nonhomologous_geometry(),
            1.0e-12,
            np.zeros((2, 2)),
            [3.95e14, 3.90e14],
            0,
            0,
            4.4139581936031586e10,
            InteractionType.ESCATTERING,
            0,
            -1,
            np.zeros((2, 2)),
            np.zeros((2, 2)),
        ),
        (
            _make_nonhomologous_geometry(r_outer_first_shell=2.0e16),
            1.0e-20,
            np.ones((2, 2)) * 100.0,
            [3.80e14, 3.70e14],
            0,
            0,
            1.9759209486510732e16,
            InteractionType.LINE,
            0,
            -1,
            np.array([[2.1375e-15, 0.0], [0.0, 0.0]]),
            np.array([[0.855, 0.0], [0.0, 0.0]]),
        ),
        (
            _make_nonhomologous_geometry(negative_velocity_gradient=True),
            1.0e-20,
            np.zeros((2, 2)),
            [3.95e14, 3.90e14],
            1,
            1,
            1.3294552658190883e14,
            InteractionType.BOUNDARY,
            2,
            1,
            np.zeros((2, 2)),
            np.zeros((2, 2)),
        ),
    ],
)
def test_nonhomologous_trace_packet_characterization(
    geometry: NumbaNonhomologousRadial1DGeometry,
    opacity_electron: float,
    tau_sobolev: np.ndarray,
    line_list_nu: list[float],
    next_line_id: int,
    prev_line_id: int,
    expected_distance: float,
    expected_interaction_type: InteractionType,
    expected_next_line_id: int,
    expected_prev_line_id: int,
    expected_j_blue: np.ndarray,
    expected_edot: np.ndarray,
) -> None:
    packet = _make_rpacket(next_line_id=next_line_id, prev_line_id=prev_line_id)
    opacity_state = _make_opacity_state(line_list_nu, tau_sobolev)
    estimators = _make_line_estimators()
    _seed_numba_random(1963)

    distance, interaction_type, delta_shell = nonhomologous_trace_packet(
        packet,
        geometry,
        opacity_state,
        estimators,
        opacity_electron,
        False,
        False,
    )

    npt.assert_allclose(distance, expected_distance)
    assert interaction_type == expected_interaction_type
    assert delta_shell == 1
    assert packet.next_line_id == expected_next_line_id
    assert packet.prev_line_id == expected_prev_line_id
    npt.assert_allclose(estimators.mean_intensity_blueward, expected_j_blue)
    npt.assert_allclose(estimators.energy_deposition_line_rate, expected_edot)


def test_nonhomologous_trace_packet_no_line_fallthrough_characterization() -> None:
    packet = _make_rpacket(next_line_id=2, prev_line_id=1)
    geometry = _make_nonhomologous_geometry(r_outer_first_shell=2.0e16)
    opacity_state = _make_opacity_state(
        [3.999e14, 3.998e14], np.zeros((2, 2))
    )
    estimators = _make_line_estimators()
    _seed_numba_random(1963)

    distance, interaction_type, delta_shell = nonhomologous_trace_packet(
        packet,
        geometry,
        opacity_state,
        estimators,
        1.0e-12,
        False,
        False,
    )

    npt.assert_allclose(distance, 4.4139581936031586e10)
    assert interaction_type == InteractionType.ESCATTERING
    assert delta_shell == 1
    assert packet.next_line_id == 2
    assert packet.prev_line_id == 1
    npt.assert_allclose(estimators.mean_intensity_blueward, np.zeros((2, 2)))
