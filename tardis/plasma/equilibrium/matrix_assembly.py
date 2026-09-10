"""Shared helpers for assembling rate matrices."""

from collections.abc import Iterable, Sequence

import numpy as np
import numpy.typing as npt
import pandas as pd


def sum_duplicate_rates(rates: pd.DataFrame) -> pd.DataFrame:
    """Sum duplicate ion transitions in a rate frame."""
    return rates.groupby(
        level=(
            "atomic_number",
            "ion_number",
            "ion_number_source",
            "ion_number_destination",
        )
    ).sum()


def normalize_rate_matrices(
    matrices: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Add balance diagonals and the population-normalization row."""
    for matrix in matrices:
        np.fill_diagonal(matrix, -np.sum(matrix, axis=0))
    matrices[:, 1, :] = 1.0
    return matrices


def sum_rate_frames(
    rate_frames: Iterable[pd.DataFrame],
    multipliers: Sequence[float] | None = None,
) -> pd.DataFrame:
    """Align rate frames on their union of indexes and add them.

    Parameters
    ----------
    rate_frames : iterable of pandas.DataFrame
        Rate frames with matching columns.
    multipliers : sequence of float, optional
        Factors applied to the corresponding rate frames before summation.

    Returns
    -------
    pandas.DataFrame
        The aligned sum of the rate frames.
    """
    frames = list(rate_frames)
    if not frames:
        raise ValueError("at least one rate frame is required")
    if multipliers is None:
        multipliers = [1.0] * len(frames)
    if len(multipliers) != len(frames):
        raise ValueError("multipliers must match the number of rate frames")

    if isinstance(frames[0].index, pd.MultiIndex):
        index = pd.MultiIndex.from_tuples(
            sorted(set().union(*(frame.index for frame in frames))),
            names=frames[0].index.names,
        )
    else:
        index = frames[0].index
        for frame in frames[1:]:
            index = index.union(frame.index, sort=False)
    return sum(
        frame.reindex(index, fill_value=0) * multiplier
        for frame, multiplier in zip(frames, multipliers, strict=True)
    )


def construct_rate_matrices(
    rates: pd.DataFrame,
    shape: tuple[int, int],
    source_level: str,
    destination_level: str,
) -> npt.NDArray[np.float64]:
    """Construct dense destination-row/source-column matrices.

    Parameters
    ----------
    rates : pandas.DataFrame
        Rates indexed by source and destination levels, with one column per
        shell.
    shape : tuple of int
        Matrix row and column dimensions.
    source_level, destination_level : str
        Index level names identifying matrix columns and rows.

    Returns
    -------
    numpy.typing.NDArray
        One dense matrix per rate-frame column.
    """
    matrices = np.zeros((len(rates.columns), *shape), dtype=float)
    source = rates.index.get_level_values(source_level).to_numpy(dtype=int)
    destination = rates.index.get_level_values(destination_level).to_numpy(
        dtype=int
    )
    for shell_idx, values in enumerate(rates.to_numpy().T):
        np.add.at(matrices[shell_idx], (destination, source), values)
    return matrices
