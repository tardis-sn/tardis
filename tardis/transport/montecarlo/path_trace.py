"""Shared helpers for radiative packet trace event candidates."""

import numpy as np
from numba import njit

from tardis.transport.montecarlo import njit_dict_no_parallel
from tardis.transport.montecarlo.transport_events import TransportEvent


@njit(**njit_dict_no_parallel)
def sample_event_optical_depth() -> float:
    """
    Sample the optical depth to the next Monte Carlo event.

    Returns
    -------
    float
        Sampled event optical depth.
    """
    return -np.log(np.random.random())


@njit(**njit_dict_no_parallel)
def select_nearest_trace_candidate(
    boundary_distance: float,
    boundary_event: int,
    line_distance: float,
    line_event: int,
    continuum_distance: float,
    continuum_event: int,
) -> tuple[int, float]:
    """
    Select the nearest precomputed trace event candidate.

    Parameters
    ----------
    boundary_distance : float
        Distance to the shell-boundary candidate.
    boundary_event : int
        Event code for the shell-boundary candidate.
    line_distance : float
        Distance to the line candidate.
    line_event : int
        Event code for the line candidate.
    continuum_distance : float
        Distance to the continuum or electron-scattering candidate.
    continuum_event : int
        Event code for the continuum or electron-scattering candidate.

    Returns
    -------
    tuple[int, float]
        Selected event code and distance.
    """
    event = boundary_event
    distance = boundary_distance

    if line_distance < distance:
        event = line_event
        distance = line_distance

    if continuum_distance < distance:
        event = continuum_event
        distance = continuum_distance

    return event, distance


@njit(**njit_dict_no_parallel)
def electron_event_candidate(
    tau_event: float, opacity_electron: float
) -> tuple[int, float]:
    """
    Build an electron-scattering event candidate.

    Parameters
    ----------
    tau_event : float
        Sampled event optical depth.
    opacity_electron : float
        Electron-scattering opacity.

    Returns
    -------
    tuple[int, float]
        Event code and distance.
    """
    return TransportEvent.ESCATTERING, tau_event / opacity_electron

