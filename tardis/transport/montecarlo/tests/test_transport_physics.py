import numba
import numpy as np
import numpy.testing as npt
import pytest

from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.model.geometry.radial1d_homologous import (
    NumbaHomologousRadial1DGeometry,
)
from tardis.transport.montecarlo.configuration.constants import C_SPEED_OF_LIGHT
from tardis.transport.montecarlo.modes.classic.packet_propagation import (
    packet_propagation as propagate_packet_homologous,
)
from tardis.transport.montecarlo.modes.iip.montecarlo_transport import (
    montecarlo_transport as transport_packets_continuum,
)
from tardis.transport.montecarlo.modes.iip.packet_propagation import (
    packet_propagation as propagate_packet_continuum,
)
from tardis.transport.montecarlo.modes.montecarlo_transport import (
    montecarlo_transport_with_vpackets as transport_packets_with_vpackets,
)
from tardis.transport.montecarlo.modes.nonhomologous.packet_propagation import (
    packet_propagation as propagate_packet_piecewise_linear,
)
from tardis.transport.montecarlo.packets.movement import (
    initialize_packet_frame,
)
from tardis.transport.montecarlo.packets.radiative_packet import RPacket
from tardis.transport.montecarlo.transport_physics import (
    resolve_transport_physics,
)


@pytest.mark.parametrize("enable_full_relativity", [False, True])
def test_resolve_homologous_transport_physics(
    homologous_radial_1d_geometry: NumbaHomologousRadial1DGeometry,
    enable_full_relativity: bool,
) -> None:
    transport_physics = resolve_transport_physics(
        homologous_radial_1d_geometry,
        enable_full_relativity,
        continuum_process_enabled=False,
    )

    assert transport_physics.propagate_packet is propagate_packet_homologous
    assert (
        transport_physics.transport_packets is transport_packets_with_vpackets
    )
    assert transport_physics.enable_full_relativity is enable_full_relativity


def test_resolve_continuum_transport_physics(
    homologous_radial_1d_geometry: NumbaHomologousRadial1DGeometry,
) -> None:
    transport_physics = resolve_transport_physics(
        homologous_radial_1d_geometry,
        enable_full_relativity=True,
        continuum_process_enabled=True,
    )

    assert transport_physics.propagate_packet is propagate_packet_continuum
    assert transport_physics.transport_packets is transport_packets_continuum
    assert transport_physics.enable_full_relativity is True


def test_partial_relativity_continuum_transport_is_rejected(
    homologous_radial_1d_geometry: NumbaHomologousRadial1DGeometry,
) -> None:
    with pytest.raises(
        NotImplementedError,
        match="Continuum transport currently requires full relativity",
    ):
        resolve_transport_physics(
            homologous_radial_1d_geometry,
            enable_full_relativity=False,
            continuum_process_enabled=True,
        )


def test_resolve_piecewise_linear_transport_physics(
    radial_1d_geometry: NumbaRadial1DGeometry,
) -> None:
    transport_physics = resolve_transport_physics(
        radial_1d_geometry,
        False,
        continuum_process_enabled=False,
    )

    assert (
        transport_physics.propagate_packet is propagate_packet_piecewise_linear
    )
    assert (
        transport_physics.transport_packets is transport_packets_with_vpackets
    )
    assert transport_physics.enable_full_relativity is False


def test_piecewise_linear_full_relativity_is_rejected(
    radial_1d_geometry: NumbaRadial1DGeometry,
) -> None:
    with pytest.raises(
        NotImplementedError,
        match=r"Full relativity.*piecewise-linear kinematics",
    ):
        resolve_transport_physics(
            radial_1d_geometry,
            enable_full_relativity=True,
            continuum_process_enabled=False,
        )


def test_piecewise_linear_continuum_transport_is_rejected(
    radial_1d_geometry: NumbaRadial1DGeometry,
) -> None:
    with pytest.raises(
        NotImplementedError,
        match=r"Continuum transport.*piecewise-linear kinematics",
    ):
        resolve_transport_physics(
            radial_1d_geometry,
            enable_full_relativity=False,
            continuum_process_enabled=True,
        )


def test_unsupported_geometry_is_rejected() -> None:
    with pytest.raises(
        TypeError, match="Unsupported transport geometry: object"
    ):
        resolve_transport_physics(
            object(),
            enable_full_relativity=False,
            continuum_process_enabled=False,
        )


@pytest.mark.parametrize("enable_full_relativity", [False, True])
def test_initialize_packet_frame_homologous(
    enable_full_relativity: bool,
) -> None:
    time_explosion = 1.0e6
    beta = 0.1
    packet_radius = beta * C_SPEED_OF_LIGHT * time_explosion
    geometry = NumbaHomologousRadial1DGeometry(
        np.array([0.05 * C_SPEED_OF_LIGHT * time_explosion]),
        np.array([0.15 * C_SPEED_OF_LIGHT * time_explosion]),
        np.array([0.05 * C_SPEED_OF_LIGHT]),
        np.array([0.15 * C_SPEED_OF_LIGHT]),
        time_explosion,
    )
    initial_mu = 0.6
    initial_nu = 1.0e15
    initial_energy = 1.0
    packet = RPacket(
        r=packet_radius,
        mu=initial_mu,
        nu=initial_nu,
        energy=initial_energy,
        seed=1963,
    )

    initialize_packet_frame(packet, geometry, enable_full_relativity)

    if enable_full_relativity:
        inverse_doppler_factor = (1.0 + initial_mu * beta) / np.sqrt(
            1.0 - beta**2
        )
        expected_mu = (initial_mu + beta) / (1.0 + beta * initial_mu)
    else:
        inverse_doppler_factor = 1.0 / (1.0 - initial_mu * beta)
        expected_mu = initial_mu
    npt.assert_allclose(packet.nu, initial_nu * inverse_doppler_factor)
    npt.assert_allclose(packet.energy, initial_energy * inverse_doppler_factor)
    npt.assert_allclose(packet.mu, expected_mu)
    if not numba.config.DISABLE_JIT:
        assert initialize_packet_frame.nopython_signatures


def test_initialize_packet_frame_uses_piecewise_linear_velocity() -> None:
    geometry = NumbaRadial1DGeometry(
        np.array([1.0e14]),
        np.array([2.0e14]),
        np.array([1.0e8]),
        np.array([2.0e8]),
    )
    initial_mu = 0.3
    initial_nu = 1.0e15
    initial_energy = 1.0
    packet = RPacket(
        r=1.25e14,
        mu=initial_mu,
        nu=initial_nu,
        energy=initial_energy,
        seed=1963,
    )
    local_velocity = 1.25e8
    inverse_doppler_factor = 1.0 / (
        1.0 - initial_mu * local_velocity / C_SPEED_OF_LIGHT
    )

    initialize_packet_frame(packet, geometry, False)

    npt.assert_allclose(packet.nu, initial_nu * inverse_doppler_factor)
    npt.assert_allclose(packet.energy, initial_energy * inverse_doppler_factor)
    assert packet.mu == initial_mu
