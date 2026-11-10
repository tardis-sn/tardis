"""Orchestration, import, and compilation contracts for explicit physics."""

import ast
import importlib
import inspect
import subprocess
import sys
from dataclasses import MISSING, FrozenInstanceError, fields, replace
from enum import IntEnum, StrEnum
from itertools import product
from pathlib import Path
from types import SimpleNamespace

import numba as nb
import numpy as np
import pytest
from astropy import units as u

from tardis.transport.montecarlo import physics as physics_module
from tardis.transport.montecarlo.configuration import base as configuration
from tardis.transport.montecarlo.kernels import configuration_enums
from tardis.transport.montecarlo.kernels.configuration_enums import (
    LineInteractionType,
    RelativityType,
)
from tardis.transport.montecarlo.physics import (
    SUPPORTED_TRANSPORT_PHYSICS,
    Kinematics,
    LineInteraction,
    NumbaTransportPhysics,
    Relativity,
    ResolvedTransportKernel,
    TransportPhysics,
    validate_transport_physics,
)

# Independent capability specification: (kinematics, relativity, lines, continuum).
EXPECTED_REGIMES = {
    ("homologous", "partial", "none", False),
    ("homologous", "partial", "scatter", False),
    ("homologous", "partial", "downbranch", False),
    ("homologous", "partial", "macroatom", False),
    ("homologous", "full", "none", False),
    ("homologous", "full", "scatter", False),
    ("homologous", "full", "downbranch", False),
    ("homologous", "full", "macroatom", False),
    ("homologous", "full", "macroatom", True),
    ("nonhomologous", "partial", "none", False),
    ("nonhomologous", "partial", "scatter", False),
    ("nonhomologous", "partial", "downbranch", False),
    ("nonhomologous", "partial", "macroatom", False),
}
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
def test_validate_all_transport_physics(physics: TransportPhysics) -> None:
    regime = (
        physics.kinematics.value,
        physics.relativity.value,
        physics.line_interaction.value,
        physics.continuum_processes,
    )
    if regime in EXPECTED_REGIMES:
        assert validate_transport_physics(physics) is None
    else:
        with pytest.raises(NotImplementedError):
            validate_transport_physics(physics)


def test_capability_set_has_exactly_26_supported_and_38_rejected() -> None:
    expected = {
        TransportPhysics(
            Kinematics(kinematics),
            Relativity(relativity),
            LineInteraction(lines),
            electron_scattering,
            continuum,
        )
        for kinematics, relativity, lines, continuum in EXPECTED_REGIMES
        for electron_scattering in (False, True)
    }
    assert isinstance(SUPPORTED_TRANSPORT_PHYSICS, frozenset)
    assert SUPPORTED_TRANSPORT_PHYSICS == expected
    assert len(expected) == 26
    assert len(set(ALL_PHYSICS) - expected) == 38


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
    physics = replace(explicit_physics, **{component: malformed})
    with pytest.raises(TypeError, match=component):
        validate_transport_physics(physics)


@pytest.mark.parametrize("malformed", [None, {}, "homologous", True])
def test_validate_requires_transport_physics(malformed: object) -> None:
    with pytest.raises(TypeError, match="TransportPhysics"):
        validate_transport_physics(malformed)


def test_transport_physics_requires_all_five_components(
    explicit_physics: TransportPhysics,
) -> None:
    physics_fields = fields(TransportPhysics)
    assert len(physics_fields) == 5
    components = {
        field.name: getattr(explicit_physics, field.name)
        for field in physics_fields
    }
    for field in physics_fields:
        assert field.default is MISSING
        assert field.default_factory is MISSING
        with pytest.raises(TypeError, match=field.name):
            TransportPhysics(
                **{
                    key: value
                    for key, value in components.items()
                    if key != field.name
                }
            )
    assert not hasattr(explicit_physics, "__dict__")
    assert not hasattr(explicit_physics, "to_numba")
    with pytest.raises(FrozenInstanceError):
        explicit_physics.relativity = Relativity.FULL


def test_host_and_compiled_selection_contract() -> None:
    for enum, values in (
        (Kinematics, {"homologous", "nonhomologous"}),
        (Relativity, {"partial", "full"}),
        (LineInteraction, {"scatter", "downbranch", "macroatom", "none"}),
    ):
        assert issubclass(enum, StrEnum)
        assert {member.value for member in enum} == values
        assert len(enum.__members__) == len(values)
    for enum, expected in (
        (RelativityType, {"PARTIAL": 0, "FULL": 1}),
        (
            LineInteractionType,
            {"SCATTER": 0, "DOWNBRANCH": 1, "MACROATOM": 2, "NONE": 3},
        ),
    ):
        assert issubclass(enum, IntEnum)
        assert {
            name: member.value for name, member in enum.__members__.items()
        } == expected
        assert len(enum) == len(expected)
    assert {
        name for name in vars(configuration_enums) if not name.startswith("_")
    } == {
        "IntEnum",
        "RelativityType",
        "LineInteractionType",
    }


@pytest.mark.parametrize(
    "relativity,line_interaction", product(RelativityType, LineInteractionType)
)
def test_numba_physics_retains_selections(
    relativity: RelativityType, line_interaction: LineInteractionType
) -> None:
    physics = NumbaTransportPhysics(relativity, line_interaction, False, True)
    assert physics.relativity is relativity
    assert physics.line_interaction is line_interaction
    assert physics.electron_scattering_enabled is False
    assert physics.continuum_processes_enabled is True


@pytest.mark.parametrize(
    "module_name",
    [
        "tardis.transport.montecarlo.interaction_events",
        "tardis.transport.montecarlo.interaction_event_callers",
        "tardis.transport.montecarlo.configuration.base",
        "tardis.transport.montecarlo.modes.nonhomologous.interaction_events",
        "tardis.transport.montecarlo.modes.nonhomologous.interaction_event_callers",
        "tardis.transport.montecarlo.modes.iip.interaction_event_callers",
    ],
)
def test_line_selection_has_one_shared_identity(module_name: str) -> None:
    assert (
        importlib.import_module(module_name).LineInteractionType
        is LineInteractionType
    )


@pytest.mark.parametrize(
    "line_interaction", ["scatter", "downbranch", "macroatom"]
)
def test_configuration_preserves_line_selection(
    line_interaction: str,
    montecarlo_configuration: configuration.MonteCarloConfiguration,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
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


def test_numba_physics_construction_and_access() -> None:
    @nb.njit
    def construct_and_read() -> tuple[
        RelativityType, LineInteractionType, bool, bool
    ]:
        physics = NumbaTransportPhysics(
            RelativityType.FULL, LineInteractionType.MACROATOM, False, True
        )
        return (
            physics.relativity,
            physics.line_interaction,
            physics.electron_scattering_enabled,
            physics.continuum_processes_enabled,
        )

    assert construct_and_read() == (
        RelativityType.FULL,
        LineInteractionType.MACROATOM,
        False,
        True,
    )

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
    if nb.config.DISABLE_JIT:
        assert not hasattr(NumbaTransportPhysics, "class_type")
    else:
        assert construct_and_read.nopython_signatures
        assert count_selected.nopython_signatures
        assert count_selected.targetoptions["parallel"]
        compiled_loop = count_selected.overloads[count_selected.signatures[0]]
        assert compiled_loop.metadata["parfor_diagnostics"].count_parfors() == 1
        assert NumbaTransportPhysics.class_type.struct == {
            "relativity": nb.typeof(RelativityType.PARTIAL),
            "line_interaction": nb.typeof(LineInteractionType.SCATTER),
            "electron_scattering_enabled": nb.boolean,
            "continuum_processes_enabled": nb.boolean,
        }
        context = nb.core.registry.cpu_target.target_context
        for enum_type in (RelativityType, LineInteractionType):
            assert context.get_data_type(
                nb.typeof(next(iter(enum_type)))
            ) == context.get_data_type(nb.int64)


def test_resolved_transport_kernel_contract() -> None:
    @nb.njit
    def dispatcher() -> None:
        return

    compiled_physics = NumbaTransportPhysics(
        RelativityType.PARTIAL, LineInteractionType.SCATTER, True, False
    )
    resolved = ResolvedTransportKernel(dispatcher, dispatcher, compiled_physics)
    assert [field.name for field in fields(resolved)] == [
        "driver",
        "packet_propagator",
        "numba_physics",
    ]
    assert resolved.driver is dispatcher
    assert resolved.packet_propagator is dispatcher
    assert resolved.numba_physics is compiled_physics
    assert not hasattr(resolved, "__dict__")
    with pytest.raises(FrozenInstanceError):
        resolved.driver = dispatcher


@pytest.mark.parametrize(
    "module",
    [
        "tardis.transport.montecarlo.physics",
        "tardis.opacities.opacity_state",
        "tardis.io.model.parse_packet_source_configuration",
    ],
)
def test_fresh_process_imports(module: str) -> None:
    command = (
        "import importlib, sys; "
        f"importlib.import_module({module!r}); "
        "assert not any(name in sys.modules for name in "
        "('tardis.transport.montecarlo.kernels.driver', "
        "'tardis.transport.montecarlo.kernels.processes', "
        "'tardis.transport.montecarlo.kernels.kinematics'))"
    )
    subprocess.run([sys.executable, "-c", command], check=True, timeout=60)


def test_physics_imports_only_leaf_configuration_enums() -> None:
    tree = ast.parse(inspect.getsource(physics_module))
    kernel_imports = [
        node.module
        for node in tree.body
        if isinstance(node, ast.ImportFrom)
        and node.module
        and ".kernels" in node.module
    ]
    assert kernel_imports == [
        "tardis.transport.montecarlo.kernels.configuration_enums"
    ]
    selection_path = Path(configuration_enums.__file__)
    imports = [
        node
        for node in ast.walk(ast.parse(inspect.getsource(configuration_enums)))
        if isinstance(node, (ast.Import, ast.ImportFrom))
    ]
    assert len(imports) == 1
    assert isinstance(imports[0], ast.ImportFrom)
    assert imports[0].module == "enum"
    initializer = selection_path.with_name("__init__.py").read_text()
    assert not any(
        isinstance(node, (ast.Import, ast.ImportFrom))
        for node in ast.walk(ast.parse(initializer))
    )
