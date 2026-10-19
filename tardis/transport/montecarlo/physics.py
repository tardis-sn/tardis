"""Explicit host physics selections and the scalar compiled representation."""

from dataclasses import dataclass
from enum import StrEnum

import numba as nb
from numba.core.registry import CPUDispatcher
from numba.experimental import jitclass

from tardis.transport.montecarlo.kernels.configuration_enums import (
    LineInteractionType,
    RelativityType,
)


class Kinematics(StrEnum):
    """Select the explicit ejecta velocity law."""

    HOMOLOGOUS = "homologous"
    NONHOMOLOGOUS = "nonhomologous"


class Relativity(StrEnum):
    """Select the relativistic transport treatment."""

    PARTIAL = "partial"
    FULL = "full"


class LineInteraction(StrEnum):
    """Select the bound-bound interaction treatment."""

    NONE = "none"
    SCATTER = "scatter"
    DOWNBRANCH = "downbranch"
    MACROATOM = "macroatom"


@jitclass(
    [
        ("relativity", nb.typeof(RelativityType.PARTIAL)),
        ("line_interaction", nb.typeof(LineInteractionType.SCATTER)),
        ("electron_scattering_enabled", nb.boolean),
        ("continuum_processes_enabled", nb.boolean),
    ]
)
class NumbaTransportPhysics:
    """Store the four scalar selections consumed by compiled transport.

    Parameters
    ----------
    relativity : RelativityType
        Compiled relativity selection.
    line_interaction : LineInteractionType
        Compiled line-interaction selection.
    electron_scattering_enabled : bool
        Whether electron scattering contributes opacity and events.
    continuum_processes_enabled : bool
        Whether continuum processes contribute opacity and events.

    Notes
    -----
    Kinematics is represented by the selected packet propagator, not a scalar.
    Numba lowers each enum field to an integer scalar while retaining its type.
    """

    def __init__(
        self,
        relativity: RelativityType,
        line_interaction: LineInteractionType,
        electron_scattering_enabled: bool,
        continuum_processes_enabled: bool,
    ) -> None:
        self.relativity = relativity
        self.line_interaction = line_interaction
        self.electron_scattering_enabled = electron_scattering_enabled
        self.continuum_processes_enabled = continuum_processes_enabled


@dataclass(frozen=True, slots=True)
class TransportPhysics:
    """Store all five explicitly required host physics selections.

    Parameters
    ----------
    kinematics : Kinematics
        Ejecta velocity law.
    relativity : Relativity
        Relativistic transport treatment.
    line_interaction : LineInteraction
        Bound-bound interaction treatment.
    electron_scattering : bool
        Whether electron scattering is enabled.
    continuum_processes : bool
        Whether continuum processes are enabled.
    """

    kinematics: Kinematics
    relativity: Relativity
    line_interaction: LineInteraction
    electron_scattering: bool
    continuum_processes: bool


SUPPORTED_TRANSPORT_PHYSICS: frozenset[TransportPhysics] = frozenset(
    {
        TransportPhysics(
            Kinematics.HOMOLOGOUS, relativity, lines, electrons, False
        )
        for relativity in Relativity
        for lines in LineInteraction
        for electrons in (False, True)
    }
    | {
        TransportPhysics(
            Kinematics.HOMOLOGOUS,
            Relativity.FULL,
            LineInteraction.MACROATOM,
            electrons,
            True,
        )
        for electrons in (False, True)
    }
    | {
        TransportPhysics(
            Kinematics.NONHOMOLOGOUS,
            Relativity.PARTIAL,
            lines,
            electrons,
            False,
        )
        for lines in LineInteraction
        for electrons in (False, True)
    }
)


def validate_transport_physics(physics: TransportPhysics) -> None:
    """Check component types and the supported real-transport capabilities.

    Parameters
    ----------
    physics : TransportPhysics
        Explicit transport selections to validate.

    Raises
    ------
    TypeError
        If the object, an enum component, or a boolean component is malformed.
    NotImplementedError
        If a well-typed combination is unsupported.
    """
    if not isinstance(physics, TransportPhysics):
        raise TypeError("physics must be a TransportPhysics")
    for name, enum_type in (
        ("kinematics", Kinematics),
        ("relativity", Relativity),
        ("line_interaction", LineInteraction),
    ):
        if not isinstance(getattr(physics, name), enum_type):
            raise TypeError(f"{name} must be a {enum_type.__name__}")
    for name in ("electron_scattering", "continuum_processes"):
        if type(getattr(physics, name)) is not bool:
            raise TypeError(f"{name} must be a bool")
    if physics not in SUPPORTED_TRANSPORT_PHYSICS:
        raise NotImplementedError(
            f"Unsupported transport physics: {physics}. "
            "Nonhomologous transport requires partial relativity and no "
            "continuum; continuum requires homologous full-relativity "
            "macroatom transport."
        )


@dataclass(frozen=True, slots=True)
class ResolvedTransportKernel:
    """Bundle selected compiled dispatchers and physics on the host only.

    Parameters
    ----------
    driver : numba.core.registry.CPUDispatcher
        Compiled real-packet transport driver.
    packet_propagator : numba.core.registry.CPUDispatcher
        Compiled propagator selected for the explicit kinematics.
    numba_physics : NumbaTransportPhysics
        Four scalar physics fields passed to compiled transport.
    """

    driver: CPUDispatcher
    packet_propagator: CPUDispatcher
    numba_physics: NumbaTransportPhysics
