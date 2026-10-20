"""Stable identity for an explorer input state."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, is_dataclass

LENGTH_UNIT = "illustrative length unit"
VELOCITY_UNIT = "illustrative velocity unit"


def input_fingerprint(inputs: object) -> str:
    """Return a stable, short identifier for exact canonical input values.

    Parameters
    ----------
    inputs : object
        Immutable explorer input dataclass.

    Returns
    -------
    str
        First 16 hexadecimal digits of a SHA-256 digest.
    """
    if not is_dataclass(inputs):
        raise TypeError("inputs must be a dataclass")
    payload = {
        "inputs": asdict(inputs),
        "units": {"length": LENGTH_UNIT, "velocity": VELOCITY_UNIT},
    }
    canonical = json.dumps(
        payload, allow_nan=False, sort_keys=True, separators=(",", ":")
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:16]
