"""Shared transport event identifiers and trace-result tuple fields."""

from enum import IntEnum

from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
)


class TransportEvent(IntEnum):
    """
    Shared radiative transport event identifiers.

    Values intentionally mirror :class:`InteractionType` while transport
    callers are migrated from interaction-specific naming.
    """

    NO_INTERACTION = InteractionType.NO_INTERACTION
    BOUNDARY = InteractionType.BOUNDARY
    LINE = InteractionType.LINE
    ESCATTERING = InteractionType.ESCATTERING
    CONTINUUM_PROCESS = InteractionType.CONTINUUM_PROCESS


TRACE_RESULT_DISTANCE = 0
TRACE_RESULT_EVENT = 1
TRACE_RESULT_DELTA_SHELL = 2

