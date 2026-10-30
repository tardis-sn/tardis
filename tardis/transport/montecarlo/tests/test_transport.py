from collections.abc import Callable

import numpy as np
import numpy.testing as npt
import pytest

from tardis import constants as const
from tardis.conftest import sync_ndarray_assert_allclose
from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.model.geometry.radial1d_homologous import (
    NumbaHomologousRadial1DGeometry,
)
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.transport.montecarlo.estimators.estimators_bulk import (
    EstimatorsBulk,
    init_estimators_bulk,
)
from tardis.transport.montecarlo.estimators.estimators_line import (
    EstimatorsLine,
)
from tardis.transport.montecarlo.modes.homologous_rad_packet_transport import (
    trace_packet as homologous_trace_packet,
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
    RPacket,
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
    homologous_radial_1d_geometry: NumbaHomologousRadial1DGeometry,
    enable_full_relativity: bool,
    regression_data,
) -> None:
    packet = parametrized_packet
    packet.current_shell_id = 1
    estimators = init_estimators_bulk(3)

    move_r_packet(
        packet,
        1.0e13,
        homologous_radial_1d_geometry,
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
    radial_1d_geometry: NumbaRadial1DGeometry,
    bulk_estimators,
    regression_data,
) -> None:
    packet = parametrized_packet
    packet.current_shell_id = 0

    move_r_packet(
        packet,
        1.0e13,
        radial_1d_geometry,
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
    homologous_radial_1d_geometry: NumbaHomologousRadial1DGeometry,
    regression_data,
) -> None:
    packet = parametrized_packet
    packet.current_shell_id = 1
    estimators = init_estimators_bulk(3)

    move_r_packet(
        packet,
        0.0,
        homologous_radial_1d_geometry,
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
    radial_1d_geometry: NumbaRadial1DGeometry,
    bulk_estimators,
    regression_data,
) -> None:
    packet = parametrized_packet
    packet.current_shell_id = 0

    move_r_packet(
        packet,
        0.0,
        radial_1d_geometry,
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
        "homologous_radial_1d_geometry",
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
    indirect=["opacity_state_args", "homologous_radial_1d_geometry"],
)
def test_classic_trace_packet(
    parametrized_packet,
    opacity_electron: float,
    classic_opacity_state,
    disable_line_scattering: bool,
    homologous_radial_1d_geometry: NumbaHomologousRadial1DGeometry,
    line_estimators,
    set_seed_fixture,
    expected_interaction_type: InteractionType,
    regression_data,
) -> None:
    set_seed_fixture(1963)

    distance, interaction_type, delta_shell = homologous_trace_packet(
        parametrized_packet,
        homologous_radial_1d_geometry,
        5.2e7,
        classic_opacity_state,
        line_estimators,
        opacity_electron,
        1.0,
        False,
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
@pytest.mark.parametrize("homologous_radial_1d_geometry", [2.0e16], indirect=True)
@pytest.mark.parametrize(
    "opacity_state_args",
    [FALLTHROUGH_OPACITY],
    indirect=True,
)
def test_classic_trace_packet_no_line_fallthrough(
    parametrized_packet,
    homologous_radial_1d_geometry: NumbaHomologousRadial1DGeometry,
    classic_opacity_state,
    line_estimators,
    set_seed_fixture,
    regression_data,
) -> None:
    # Start beyond the available line list to lock in the no-line fallthrough.
    set_seed_fixture(1963)

    distance, interaction_type, delta_shell = homologous_trace_packet(
        parametrized_packet,
        homologous_radial_1d_geometry,
        5.2e7,
        classic_opacity_state,
        line_estimators,
        1.0e-12,
        1.0,
        False,
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
        "homologous_radial_1d_geometry",
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
    indirect=["opacity_state_args", "homologous_radial_1d_geometry"],
)
def test_iip_trace_packet(
    parametrized_packet,
    chi_continuum: float,
    escat_prob: float,
    iip_opacity_state,
    disable_line_scattering: bool,
    homologous_radial_1d_geometry: NumbaHomologousRadial1DGeometry,
    line_estimators,
    set_seed_fixture,
    expected_interaction_type: InteractionType,
    regression_data,
) -> None:
    set_seed_fixture(1963)

    distance, interaction_type, delta_shell = homologous_trace_packet(
        parametrized_packet,
        homologous_radial_1d_geometry,
        5.2e7,
        iip_opacity_state,
        line_estimators,
        chi_continuum,
        escat_prob,
        True,
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
@pytest.mark.parametrize("homologous_radial_1d_geometry", [2.0e16], indirect=True)
@pytest.mark.parametrize(
    "opacity_state_args",
    [FALLTHROUGH_OPACITY],
    indirect=True,
)
def test_iip_trace_packet_no_line_fallthrough(
    parametrized_packet,
    homologous_radial_1d_geometry: NumbaHomologousRadial1DGeometry,
    iip_opacity_state,
    line_estimators,
    set_seed_fixture,
    regression_data,
) -> None:
    # Start beyond the available line list to lock in the no-line fallthrough.
    set_seed_fixture(1963)

    distance, interaction_type, delta_shell = homologous_trace_packet(
        parametrized_packet,
        homologous_radial_1d_geometry,
        5.2e7,
        iip_opacity_state,
        line_estimators,
        1.0e-12,
        1.0,
        True,
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
        "radial_1d_geometry",
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
            InteractionType.LINE,
            0,
            -1,
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
        "radial_1d_geometry",
        "opacity_state_args",
        "parametrized_packet",
    ],
)
def test_nonhomologous_trace_packet(
    radial_1d_geometry: NumbaRadial1DGeometry,
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

    distance, interaction_type, delta_shell, _, _ = nonhomologous_trace_packet(
        parametrized_packet,
        radial_1d_geometry,
        classic_opacity_state,
        line_estimators,
        opacity_electron,
        False,
        False,
        0,
        0,
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
    "radial_1d_geometry",
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
    radial_1d_geometry: NumbaRadial1DGeometry,
    classic_opacity_state,
    line_estimators,
    set_seed_fixture,
    regression_data,
) -> None:
    # Start beyond the available line list to lock in the no-line fallthrough.
    set_seed_fixture(1963)

    distance, interaction_type, delta_shell, _, _ = nonhomologous_trace_packet(
        parametrized_packet,
        radial_1d_geometry,
        classic_opacity_state,
        line_estimators,
        1.0e-12,
        False,
        False,
        0,
        0,
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


@pytest.fixture
def shell_gradient_reversal_geometry() -> NumbaRadial1DGeometry:
    """Return continuous shells with equal and opposite velocity gradients."""
    return NumbaRadial1DGeometry(
        np.array([1.0e14, 1.5e14]),
        np.array([1.5e14, 2.0e14]),
        np.array([1.0e8, 1.4e8]),
        np.array([1.4e8, 1.0e8]),
    )


@pytest.mark.parametrize(
    "line_search_packet", [{"r": 1.1e14, "mu": 1.0}], indirect=True
)
def test_nonhomologous_shell_gradient_reversal_revisits_line(
    line_search_packet: RPacket,
    shell_gradient_reversal_geometry: NumbaRadial1DGeometry,
    classic_opacity_state: OpacityStateNumba,
    bulk_estimators: EstimatorsBulk,
    line_estimators: EstimatorsLine,
    set_seed_fixture: Callable[[int], None],
) -> None:
    """A continuous shell edge can reverse drift and expose a line again.

    Claim: The same line contributes once on each side of a gradient reversal.
    Regime: Radial ray through two velocity-continuous, linear shells.
    Verification: Both shell profiles reach the same projected velocity at
    radii 1.25e14 and 1.75e14 cm.
    """
    geometry = shell_gradient_reversal_geometry
    packet = line_search_packet
    estimators = line_estimators
    line_frequency = packet.nu * (1.0 - 1.2e8 / const.c.cgs.value)
    classic_opacity_state.line_list_nu = np.array([line_frequency])
    set_seed_fixture(1963)

    (
        first_boundary_distance,
        first_event,
        first_shell_step,
        next_untraversed_line_id,
        line_search_direction,
    ) = nonhomologous_trace_packet(
        packet,
        geometry,
        classic_opacity_state,
        estimators,
        1.0e-30,
        False,
        True,
        0,
        0,
    )
    assert first_event == InteractionType.BOUNDARY
    npt.assert_allclose(first_boundary_distance, 4.0e13, rtol=0.0, atol=1.0)
    npt.assert_allclose(
        estimators.energy_deposition_line_rate[0, 0],
        line_frequency / packet.nu,
        rtol=1.0e-12,
    )

    move_r_packet(packet, first_boundary_distance, geometry, bulk_estimators, False)
    move_packet_across_shell_boundary(packet, first_shell_step, 2)
    assert packet.current_shell_id == 1

    second_boundary_distance, second_event, _, _, _ = nonhomologous_trace_packet(
        packet,
        geometry,
        classic_opacity_state,
        estimators,
        1.0e-30,
        False,
        True,
        next_untraversed_line_id,
        line_search_direction,
    )
    assert second_event == InteractionType.BOUNDARY
    npt.assert_allclose(second_boundary_distance, 5.0e13, rtol=0.0, atol=1.0)
    npt.assert_allclose(
        estimators.energy_deposition_line_rate[0, 1],
        line_frequency / packet.nu,
        rtol=1.0e-12,
    )
