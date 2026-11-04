"""Prepare epoch-specific volumetric deposition from high-energy transport."""

from dataclasses import dataclass
from typing import Protocol

import numpy as np
import pandas as pd
from astropy import units as u


class HighEnergyDepositionResult(Protocol):
    """High-energy outputs required to prepare a deposition state."""

    gamma_ray_deposition_estimator: pd.DataFrame
    positron_energy: pd.DataFrame


@dataclass(frozen=True)
class DepositionState:
    """Deposition rate densities for one high-energy time bin."""

    epoch: u.Quantity
    time_bin_idx: int
    shell_index: pd.Index
    gamma_rate_density: u.Quantity
    positron_rate_density: u.Quantity
    rate_density: u.Quantity
    includes_pair_creation: bool = False


def prepare_deposition_state(
    high_energy_result: HighEnergyDepositionResult,
    time_bin_edges: u.Quantity,
    epoch: u.Quantity,
    shell_volume: u.Quantity,
    shell_index: pd.Index,
) -> DepositionState:
    """Select an epoch and combine flight-path and positron deposition.

    Parameters
    ----------
    high_energy_result : HighEnergyDepositionResult
        High-energy transport deposition outputs.
    time_bin_edges : astropy.units.Quantity
        Boundaries of the high-energy transport time bins.
    epoch : astropy.units.Quantity
        Optical epoch used to select the containing high-energy time bin.
    shell_volume : astropy.units.Quantity
        Shell volumes with shape ``(shell, time_bin)``.
    shell_index : pandas.Index
        Required optical shell labels and ordering.

    Returns
    -------
    DepositionState
        Aligned gamma, positron, and combined rate densities.
    """
    time_bin_edges = time_bin_edges.to_value(u.s)
    if np.any(np.diff(time_bin_edges) <= 0.0):
        raise ValueError("High-energy time-bin edges must increase.")
    epoch_seconds = epoch.to_value(u.s)
    if epoch_seconds < time_bin_edges[0] or epoch_seconds > time_bin_edges[-1]:
        raise ValueError("The optical epoch is outside the high-energy grid.")
    time_bin_idx = min(
        np.searchsorted(time_bin_edges, epoch_seconds, side="right") - 1,
        len(time_bin_edges) - 2,
    )

    gamma = high_energy_result.gamma_ray_deposition_estimator
    positron = high_energy_result.positron_energy
    if set(gamma.index) != set(shell_index) or set(positron.index) != set(
        shell_index
    ):
        raise ValueError(
            "High-energy and optical shell labels do not describe the same grid."
        )
    if len(gamma.columns) != len(
        time_bin_edges
    ) - 1 or not gamma.columns.equals(positron.columns):
        raise ValueError("High-energy deposition time bins do not align.")
    expected_columns = pd.Index(time_bin_edges[:-1])
    if not np.allclose(
        gamma.columns.to_numpy(dtype=float),
        expected_columns.to_numpy(dtype=float),
        rtol=0.0,
        atol=0.0,
    ):
        raise ValueError(
            "High-energy output columns do not match time-bin edges."
        )

    shell_volume = shell_volume.to(u.cm**3)
    expected_shape = (len(shell_index), len(time_bin_edges) - 1)
    if shell_volume.shape != expected_shape:
        raise ValueError(
            "shell_volume must have shape (number of shells, number of time bins)."
        )
    bin_duration = (
        time_bin_edges[time_bin_idx + 1] - time_bin_edges[time_bin_idx]
    ) * u.s
    gamma_rate_density = (
        gamma.reindex(shell_index).iloc[:, time_bin_idx].to_numpy(dtype=float)
        * u.erg
        / u.s
        / u.cm**3
    )
    positron_rate_density = (
        positron.reindex(shell_index)
        .iloc[:, time_bin_idx]
        .to_numpy(dtype=float)
        * u.erg
        / bin_duration
        / shell_volume[:, time_bin_idx]
    ).to(u.erg / u.s / u.cm**3)
    if np.any(~np.isfinite(gamma_rate_density.value)) or np.any(
        gamma_rate_density.value < 0.0
    ):
        raise ValueError(
            "The flight-path gamma deposition must be finite and nonnegative."
        )
    if np.any(~np.isfinite(positron_rate_density.value)) or np.any(
        positron_rate_density.value < 0.0
    ):
        raise ValueError(
            "The positron deposition must be finite and nonnegative."
        )
    return DepositionState(
        epoch.to(u.s),
        time_bin_idx,
        shell_index.copy(),
        gamma_rate_density,
        positron_rate_density,
        gamma_rate_density + positron_rate_density,
    )
