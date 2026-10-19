"""Shared integer physics enums for compiled transport and its callers."""

from enum import IntEnum


class RelativityType(IntEnum):
    """Select the relativistic treatment in compiled transport."""

    PARTIAL = 0
    FULL = 1


class LineInteractionType(IntEnum):
    """Select bound-bound interactions with the existing integer semantics."""

    SCATTER = 0
    DOWNBRANCH = 1
    MACROATOM = 2
    NONE = 3
