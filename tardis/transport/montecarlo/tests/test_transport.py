import numpy as np
import pytest

from tardis.conftest import sync_ndarray_assert_allclose
from tardis.transport.montecarlo.estimators.estimators_bulk import (
    init_estimators_bulk,
)
from tardis.transport.montecarlo.modes.classic.rad_packet_transport import (
    trace_packet as classic_trace_packet,
)
from tardis.transport.montecarlo.modes.iip.rad_packet_transport import (
    trace_packet as iip_trace_packet,
)
from tardis.transport.montecarlo.modes.nonhomologous.rad_packet_transport import (
    trace_packet as nonhomologous_trace_packet,
)
from tardis.transport.montecarlo.packets.movement import (
    move_packet_across_shell_boundary,
    move_r_packet,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
    PacketStatus,
)

RTOL = 1.0e-12

NO_LINE_OPACITY = {
    "tau_sobolev": np.zeros((2, 2)),
    "line_list_nu": [3.95e14, 3.90e14],
}
LINE_OPACITY = {
    "tau_sobolev": np.ones((2, 2)) * 100.0,
    "line_list_nu": [3.999e14, 3.998e14],
}
NONHOMOLOGOUS_LINE_OPACITY = {
    "tau_sobolev": np.ones((2, 2)) * 100.0,
    "line_list_nu": [3.80e14, 3.70e14],
}
FALLTHROUGH_OPACITY = {
    "tau_sobolev": np.zeros((2, 2)),
    "line_list_nu": [3.999e14, 3.998e14],
}


@pytest.mark.parametrize("enable_full_relativity", [False, True])
def test_homologous_move_r_packet(
    parametrized_packet,
    radial_geometry,
    enable_full_relativity: bool,
    regression_data,
) -> None:
    packet = parametrized_packet
    packet.current_shell_id = 1
    estimators = init_estimators_bulk(3)

    move_r_packet(
        packet,
        1.0e13,
        radial_geometry,
        estimators,
        enable_full_relativity,
    )

    sync_ndarray_assert_allclose(
        regression_data,
        np.array([packet.r, packet.mu]),
        estimators.mean_intensity_total,
        estimators.mean_frequency,
        rtol=RTOL,
    )


def test_nonhomologous_move_r_packet(
    parametrized_packet,
    nonhomologous_geometry,
    bulk_estimators,
    regression_data,
) -> None:
    packet = parametrized_packet
    packet.current_shell_id = 0

    move_r_packet(
        packet,
        1.0e13,
        nonhomologous_geometry,
        bulk_estimators,
        False,
    )

    sync_ndarray_assert_allclose(
        regression_data,
        np.array([packet.r, packet.mu]),
        bulk_estimators.mean_intensity_total,
        bulk_estimators.mean_frequency,
        rtol=RTOL,
    )


def test_homologous_move_r_packet_zero_distance(
    parametrized_packet,
    radial_geometry,
    regression_data,
) -> None:
    packet = parametrized_packet
    packet.current_shell_id = 1
    estimators = init_estimators_bulk(3)

    move_r_packet(
        packet,
        0.0,
        radial_geometry,
        estimators,
        False,
    )

    sync_ndarray_assert_allclose(
        regression_data,
        np.array([packet.r, packet.mu]),
        estimators.mean_intensity_total,
        estimators.mean_frequency,
        rtol=RTOL,
    )


def test_nonhomologous_move_r_packet_zero_distance(
    parametrized_packet,
    nonhomologous_geometry,
    bulk_estimators,
    regression_data,
) -> None:
    packet = parametrized_packet
    packet.current_shell_id = 0

    move_r_packet(
        packet,
        0.0,
        nonhomologous_geometry,
        bulk_estimators,
        False,
    )

    sync_ndarray_assert_allclose(
        regression_data,
        np.array([packet.r, packet.mu]),
        bulk_estimators.mean_intensity_total,
        bulk_estimators.mean_frequency,
        rtol=RTOL,
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
def test_move_packet_across_shell_boundary(
    parametrized_packet,
    current_shell_id: int,
    delta_shell: int,
    no_of_shells: int,
    expected_status: PacketStatus,
    expected_shell_id: int,
) -> None:
    packet = parametrized_packet
    packet.current_shell_id = current_shell_id

    move_packet_across_shell_boundary(packet, delta_shell, no_of_shells)

    assert packet.status == expected_status
    assert packet.current_shell_id == expected_shell_id


@pytest.mark.parametrize(
    (
        "opacity_electron",
        "opacity_state_args",
        "disable_line_scattering",
        "radial_geometry",
        "expected_interaction_type",
    ),
    [
        (
            1.0e-20,
            NO_LINE_OPACITY,
            False,
            8.0e14,
            InteractionType.BOUNDARY,
        ),
        (
            1.0e-12,
            NO_LINE_OPACITY,
            False,
            8.0e14,
            InteractionType.ESCATTERING,
        ),
        (
            1.0e-20,
            LINE_OPACITY,
            False,
            2.0e16,
            InteractionType.LINE,
        ),
        (
            1.0e-20,
            LINE_OPACITY,
            True,
            2.0e16,
            InteractionType.ESCATTERING,
        ),
    ],
    indirect=["opacity_state_args", "radial_geometry"],
)
def test_classic_trace_packet(
    parametrized_packet,
    opacity_electron: float,
    classic_opacity_state,
    disable_line_scattering: bool,
    radial_geometry,
    line_estimators,
    set_seed_fixture,
    expected_interaction_type: InteractionType,
    regression_data,
) -> None:
    set_seed_fixture(1963)

    distance, interaction_type, delta_shell = classic_trace_packet(
        parametrized_packet,
        radial_geometry,
        5.2e7,
        classic_opacity_state,
        line_estimators,
        opacity_electron,
        False,
        disable_line_scattering,
    )

    assert interaction_type == expected_interaction_type
    assert delta_shell == 1
    assert parametrized_packet.next_line_id == (
        1 if disable_line_scattering else 0
    )
    sync_ndarray_assert_allclose(
        regression_data,
        distance,
        line_estimators.mean_intensity_blueward,
        line_estimators.energy_deposition_line_rate,
        rtol=RTOL,
    )


@pytest.mark.parametrize(
    "parametrized_packet", [{"next_line_id": 2}], indirect=True
)
@pytest.mark.parametrize("radial_geometry", [2.0e16], indirect=True)
@pytest.mark.parametrize(
    "opacity_state_args",
    [FALLTHROUGH_OPACITY],
    indirect=True,
)
def test_classic_trace_packet_no_line_fallthrough(
    parametrized_packet,
    radial_geometry,
    classic_opacity_state,
    line_estimators,
    set_seed_fixture,
    regression_data,
) -> None:
    # Start beyond the available line list to lock in the no-line fallthrough.
    set_seed_fixture(1963)

    distance, interaction_type, delta_shell = classic_trace_packet(
        parametrized_packet,
        radial_geometry,
        5.2e7,
        classic_opacity_state,
        line_estimators,
        1.0e-12,
        False,
        False,
    )

    assert interaction_type == InteractionType.ESCATTERING
    assert delta_shell == 1
    assert parametrized_packet.next_line_id == 2
    sync_ndarray_assert_allclose(
        regression_data,
        distance,
        line_estimators.mean_intensity_blueward,
        rtol=RTOL,
    )


@pytest.mark.parametrize(
    (
        "chi_continuum",
        "escat_prob",
        "opacity_state_args",
        "disable_line_scattering",
        "radial_geometry",
        "expected_interaction_type",
    ),
    [
        (
            1.0e-20,
            0.5,
            NO_LINE_OPACITY,
            False,
            8.0e14,
            InteractionType.BOUNDARY,
        ),
        (
            1.0e-12,
            1.0,
            NO_LINE_OPACITY,
            False,
            8.0e14,
            InteractionType.ESCATTERING,
        ),
        (
            1.0e-12,
            0.0,
            NO_LINE_OPACITY,
            False,
            8.0e14,
            InteractionType.CONTINUUM_PROCESS,
        ),
        (
            1.0e-20,
            0.5,
            LINE_OPACITY,
            False,
            2.0e16,
            InteractionType.LINE,
        ),
        (
            1.0e-20,
            0.5,
            LINE_OPACITY,
            True,
            2.0e16,
            InteractionType.ESCATTERING,
        ),
    ],
    indirect=["opacity_state_args", "radial_geometry"],
)
def test_iip_trace_packet(
    parametrized_packet,
    chi_continuum: float,
    escat_prob: float,
    iip_opacity_state,
    disable_line_scattering: bool,
    radial_geometry,
    line_estimators,
    set_seed_fixture,
    expected_interaction_type: InteractionType,
    regression_data,
) -> None:
    set_seed_fixture(1963)

    distance, interaction_type, delta_shell = iip_trace_packet(
        parametrized_packet,
        radial_geometry,
        5.2e7,
        iip_opacity_state,
        line_estimators,
        chi_continuum,
        escat_prob,
        False,
        disable_line_scattering,
    )

    assert interaction_type == expected_interaction_type
    assert delta_shell == 1
    assert parametrized_packet.next_line_id == (
        1 if disable_line_scattering else 0
    )
    sync_ndarray_assert_allclose(
        regression_data,
        distance,
        line_estimators.mean_intensity_blueward,
        line_estimators.energy_deposition_line_rate,
        rtol=RTOL,
    )


@pytest.mark.parametrize(
    "parametrized_packet", [{"next_line_id": 2}], indirect=True
)
@pytest.mark.parametrize("radial_geometry", [2.0e16], indirect=True)
@pytest.mark.parametrize(
    "opacity_state_args",
    [FALLTHROUGH_OPACITY],
    indirect=True,
)
def test_iip_trace_packet_no_line_fallthrough(
    parametrized_packet,
    radial_geometry,
    iip_opacity_state,
    line_estimators,
    set_seed_fixture,
    regression_data,
) -> None:
    # Start beyond the available line list to lock in the no-line fallthrough.
    set_seed_fixture(1963)

    distance, interaction_type, delta_shell = iip_trace_packet(
        parametrized_packet,
        radial_geometry,
        5.2e7,
        iip_opacity_state,
        line_estimators,
        1.0e-12,
        1.0,
        False,
        False,
    )

    assert interaction_type == InteractionType.ESCATTERING
    assert delta_shell == 1
    assert parametrized_packet.next_line_id == 2
    sync_ndarray_assert_allclose(
        regression_data,
        distance,
        line_estimators.mean_intensity_blueward,
        rtol=RTOL,
    )


@pytest.mark.parametrize(
    (
        "nonhomologous_geometry",
        "opacity_electron",
        "opacity_state_args",
        "parametrized_packet",
        "expected_interaction_type",
        "expected_next_line_id",
        "expected_prev_line_id",
    ),
    [
        (
            {"negative_velocity_gradient": False},
            1.0e-20,
            NO_LINE_OPACITY,
            {"next_line_id": 0, "prev_line_id": 0},
            InteractionType.BOUNDARY,
            1,
            0,
        ),
        (
            {},
            1.0e-12,
            NO_LINE_OPACITY,
            {"next_line_id": 0, "prev_line_id": 0},
            InteractionType.ESCATTERING,
            1,
            0,
        ),
        (
            {"r_outer_first_shell": 2.0e16},
            1.0e-20,
            NONHOMOLOGOUS_LINE_OPACITY,
            {"next_line_id": 0, "prev_line_id": 0},
            InteractionType.BOUNDARY, # end of line list
            1,
            0,
        ),
        (
            {"negative_velocity_gradient": True},
            1.0e-20,
            NO_LINE_OPACITY,
            {"next_line_id": 1, "prev_line_id": 1},
            InteractionType.BOUNDARY,
            1,
            0,
        ),
    ],
    indirect=[
        "nonhomologous_geometry",
        "opacity_state_args",
        "parametrized_packet",
    ],
)
def test_nonhomologous_trace_packet(
    nonhomologous_geometry,
    opacity_electron: float,
    parametrized_packet,
    classic_opacity_state,
    line_estimators,
    set_seed_fixture,
    expected_interaction_type: InteractionType,
    expected_next_line_id: int,
    expected_prev_line_id: int,
    regression_data,
) -> None:
    set_seed_fixture(1963)

    distance, interaction_type, delta_shell = nonhomologous_trace_packet(
        parametrized_packet,
        nonhomologous_geometry,
        classic_opacity_state,
        line_estimators,
        opacity_electron,
        False,
        False,
    )

    assert interaction_type == expected_interaction_type
    assert delta_shell == 1
    assert parametrized_packet.next_line_id == expected_next_line_id
    assert parametrized_packet.prev_line_id == expected_prev_line_id
    sync_ndarray_assert_allclose(
        regression_data,
        distance,
        line_estimators.mean_intensity_blueward,
        line_estimators.energy_deposition_line_rate,
        rtol=RTOL,
    )


@pytest.mark.parametrize(
    "parametrized_packet",
    [{"next_line_id": 2, "prev_line_id": 1}],
    indirect=True,
)
@pytest.mark.parametrize(
    "nonhomologous_geometry",
    [{"r_outer_first_shell": 2.0e16}],
    indirect=True,
)
@pytest.mark.parametrize(
    "opacity_state_args",
    [FALLTHROUGH_OPACITY],
    indirect=True,
)
def test_nonhomologous_trace_packet_no_line_fallthrough(
    parametrized_packet,
    nonhomologous_geometry,
    classic_opacity_state,
    line_estimators,
    set_seed_fixture,
    regression_data,
) -> None:
    # Start beyond the available line list to lock in the no-line fallthrough.
    set_seed_fixture(1963)

    distance, interaction_type, delta_shell = nonhomologous_trace_packet(
        parametrized_packet,
        nonhomologous_geometry,
        classic_opacity_state,
        line_estimators,
        1.0e-12,
        False,
        False,
    )

    assert interaction_type == InteractionType.ESCATTERING
    assert delta_shell == 1
    assert parametrized_packet.next_line_id == 2
    assert parametrized_packet.prev_line_id == 1
    sync_ndarray_assert_allclose(
        regression_data,
        distance,
        line_estimators.mean_intensity_blueward,
        rtol=RTOL,
    )
