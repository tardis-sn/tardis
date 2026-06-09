import numpy as np
import numpy.testing as npt
import pytest

from tardis.model.geometry.radial1d_nonhomologous import (
    NumbaNonhomologousRadial1DGeometry,
)
from tardis.transport.montecarlo import RPacket
from tardis.transport.montecarlo.estimators.estimators_bulk import (
    init_estimators_bulk,
)
from tardis.transport.montecarlo.modes.classic.rad_packet_transport import (
    move_packet_across_shell_boundary as classic_move_boundary,
)
from tardis.transport.montecarlo.modes.classic.rad_packet_transport import (
    move_r_packet as classic_move_r_packet,
)
from tardis.transport.montecarlo.modes.iip.rad_packet_transport import (
    move_packet_across_shell_boundary as iip_move_boundary,
)
from tardis.transport.montecarlo.modes.iip.rad_packet_transport import (
    move_r_packet as iip_move_r_packet,
)
from tardis.transport.montecarlo.modes.nonhomologous.rad_packet_transport import (
    move_packet_across_shell_boundary as nonhomologous_move_boundary,
)
from tardis.transport.montecarlo.modes.nonhomologous.rad_packet_transport import (
    move_r_packet as nonhomologous_move_r_packet,
)
from tardis.transport.montecarlo.packets.radiative_packet import PacketStatus


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
        match="Full relativity not implemented for non-homologous mode.",
    ):
        nonhomologous_move_r_packet(
            packet,
            1.0e13,
            geometry,
            estimators,
            True,
        )


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
