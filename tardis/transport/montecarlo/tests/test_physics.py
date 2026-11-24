"""Validation and compilation behavior for explicit transport physics."""

from dataclasses import replace
from itertools import product
from types import SimpleNamespace

import numba as nb
import numpy as np
import pytest
from astropy import units as u

from tardis.transport.montecarlo.configuration import base as configuration
from tardis.transport.montecarlo.kernels.configuration_enums import (
    LineInteractionType,
    RelativityType,
)
from tardis.transport.montecarlo.physics import (
    Kinematics,
    LineInteraction,
    NumbaTransportPhysics,
    Relativity,
    TransportPhysics,
    validate_transport_physics,
)

ALL_PHYSICS = tuple(
    TransportPhysics(*components)
    for components in product(
        Kinematics, Relativity, LineInteraction, (False, True), (False, True)
    )
)


@pytest.fixture
def explicit_physics() -> TransportPhysics:
    return TransportPhysics(
        Kinematics.HOMOLOGOUS,
        Relativity.PARTIAL,
        LineInteraction.SCATTER,
        True,
        False,
    )


@pytest.mark.parametrize("physics", ALL_PHYSICS)
def test_transport_physics_support_matches_capability_rules(
    physics: TransportPhysics,
) -> None:
    """Claim: transport accepts exactly the implemented physical regimes.

    Regime: all host selector combinations, including both scattering flags.
    Verification: independent capability rules, exhaustively enumerated.
    """
    supported = (
        physics.kinematics is Kinematics.HOMOLOGOUS
        and (
            not physics.continuum_processes
            or (
                physics.relativity is Relativity.FULL
                and physics.line_interaction is LineInteraction.MACROATOM
            )
        )
    ) or (
        physics.kinematics is Kinematics.NONHOMOLOGOUS
        and physics.relativity is Relativity.PARTIAL
        and not physics.continuum_processes
    )
    if supported:
        assert validate_transport_physics(physics) is None
    else:
        with pytest.raises(NotImplementedError):
            validate_transport_physics(physics)


@pytest.mark.parametrize(
    "component,malformed",
    [
        ("kinematics", "homologous"),
        ("kinematics", Relativity.PARTIAL),
        ("relativity", "partial"),
        ("relativity", True),
        ("line_interaction", "scatter"),
        ("line_interaction", 0),
        ("electron_scattering", 1),
        ("electron_scattering", np.bool_(True)),
        ("continuum_processes", 0),
        ("continuum_processes", np.bool_(False)),
        ("continuum_processes", []),
    ],
)
def test_validate_malformed_components(
    explicit_physics: TransportPhysics, component: str, malformed: object
) -> None:
    """Claim: each selector rejects values outside its declared domain.

    Regime: homologous partial-relativity line scattering without continuum.
    Verification: exact field-specific type failure after one-field mutation.
    """
    physics = replace(explicit_physics, **{component: malformed})
    with pytest.raises(TypeError, match=component):
        validate_transport_physics(physics)


@pytest.mark.parametrize("malformed", [None, {}, "homologous", True])
def test_validate_requires_transport_physics(malformed: object) -> None:
    """Claim: validation rejects objects without the physics contract.

    Regime: representative scalar and container inputs.
    Verification: exact top-level type failure.
    """
    with pytest.raises(TypeError, match="TransportPhysics"):
        validate_transport_physics(malformed)


@pytest.mark.parametrize(
    "line_interaction", ["scatter", "downbranch", "macroatom"]
)
def test_configuration_preserves_line_selection(
    line_interaction: str,
    montecarlo_configuration: configuration.MonteCarloConfiguration,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Claim: configuration preserves the requested line interaction.

    Regime: homologous partial-relativity transport for each configured mode.
    Verification: independent mapping to the compiled selector values.
    """
    monkeypatch.setattr(
        configuration.montecarlo_globals, "ENABLE_RPACKET_TRACKING", False
    )
    transport = SimpleNamespace(
        line_interaction_type=line_interaction,
        enable_full_relativity=False,
        packet_source=SimpleNamespace(base_seed=1963),
        vpacket_spawn_range=SimpleNamespace(start=1.0 * u.Hz, end=2.0 * u.Hz),
        enable_vpacket_tracking=False,
        enable_rpacket_tracking=False,
    )
    configuration.configuration_initialize(montecarlo_configuration, transport)
    assert (
        montecarlo_configuration.LINE_INTERACTION_TYPE
        == {"scatter": 0, "downbranch": 1, "macroatom": 2}[line_interaction]
    )


def test_numba_physics_in_compiled_parallel_loop() -> None:
    """Claim: compiled parallel code reads every selected physics field.

    Regime: full-relativity macroatom continuum without electron scattering.
    Verification: every iteration satisfies all four exact selector checks.
    """

    @nb.njit(parallel=True)
    def count_selected(
        physics: NumbaTransportPhysics, packet_count: int
    ) -> int:
        matches = 0
        for _ in nb.prange(packet_count):
            if (
                physics.relativity == RelativityType.FULL
                and physics.line_interaction == LineInteractionType.MACROATOM
                and not physics.electron_scattering_enabled
                and physics.continuum_processes_enabled
            ):
                matches += 1
        return matches

    physics = NumbaTransportPhysics(
        RelativityType.FULL, LineInteractionType.MACROATOM, False, True
    )
    assert count_selected(physics, 8) == 8
