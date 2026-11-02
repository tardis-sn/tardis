"""Comprehensive verification of element-selected macro-atom probabilities."""

import numpy as np
import pandas as pd
import pytest

from tardis.opacities.macro_atom.absorbing_markov_chain import (
    create_absorbing_probs,
)


@pytest.fixture
def pool_network() -> tuple[pd.DataFrame, pd.DataFrame]:
    """H/He/Fe levels, terminal H, shared pools, and two different cells."""
    levels = {
        0: (1, 0, 0),
        2: (1, 1, 0),
        4: ("i", -99, -99),
        5: (2, 0, 0),
        6: ("k", -99, -99),
        9: (1, 1, 1),
        12: (26, 1, 0),
    }
    sources = [0, 0, 2, 2, 2, 4, 4, 5, 5, 6, 6, 6, 9, 12, 12]
    destinations = [2, -1, 4, 9, -1, 0, -1, 6, -1, 0, 12, -1, -1, 6, -1]
    metadata = pd.DataFrame(
        {
            "source_level_idx": sources,
            "destination_level_idx": destinations,
            "transition_type": [1 if idx >= 0 else -1 for idx in destinations],
            "source": [levels[idx] for idx in sources],
            "destination": [
                levels.get(idx, ("r", -99, -99)) for idx in destinations
            ],
        }
    )
    probabilities = pd.DataFrame(
        {
            10: [
                0.0,
                1.0,
                0.3,
                0.2,
                0.5,
                0.5,
                0.5,
                0.3,
                0.7,
                0.1,
                0.1,
                0.8,
                1.0,
                0.25,
                0.75,
            ],
            30: [
                1.0,
                0.0,
                0.2,
                0.1,
                0.7,
                0.8,
                0.2,
                0.6,
                0.4,
                0.05,
                0.05,
                0.9,
                1.0,
                0.4,
                0.6,
            ],
        }
    )
    return probabilities, metadata


@pytest.mark.parametrize(
    "mode,elements",
    [
        ("full", []),
        ("partial", ["H"]),
        ("partial", ["He"]),
        ("partial", ["Fe", "He", "H", "H"]),
        ("explicit", []),
    ],
)
def test_final_emission_law(
    pool_network: tuple[pd.DataFrame, pd.DataFrame],
    mode: str,
    elements: list[str],
) -> None:
    """Selection preserves every eventual emission channel in a mixed network.

    Regime: two cells, multiple ions, gapped indices including zero, an
    emission-only destination, zero residual/internal rows, and pool re-entry.
    Verify the mask and maps directly, then compare final channel probabilities
    with a dense solve of the original network. The reference never uses the
    reduced matrix. A 2e-14 tolerance allows floating-point LU roundoff.
    """
    probabilities, metadata = pool_network
    result = create_absorbing_probs(
        probabilities, metadata, mode=mode, elements=elements, max_workers=2
    )
    serial = create_absorbing_probs(
        probabilities, metadata, mode=mode, elements=elements
    )
    np.testing.assert_array_equal(result.matrix, serial.matrix)

    states = np.unique(metadata.source_level_idx)
    global_size = int(states.max()) + 1
    atomic_numbers = {"H": 1, "He": 2, "Fe": 26}
    selected_elements = {atomic_numbers[symbol] for symbol in elements}
    selected = (
        np.arange(global_size)
        if mode == "full"
        else np.unique(
            metadata.loc[
                metadata.source.map(
                    lambda level: level[0] in selected_elements
                ),
                "source_level_idx",
            ]
        )
    )
    lookup = np.full(global_size, -1, dtype=np.int64)
    lookup[selected] = np.arange(len(selected))
    np.testing.assert_array_equal(result.reduced_to_macro_level_idx, selected)
    np.testing.assert_array_equal(result.macro_level_to_reduced_idx, lookup)
    eliminated = (
        (metadata.transition_type >= 0)
        & metadata.source_level_idx.isin(selected)
        & metadata.destination_level_idx.isin(selected)
    )
    expected_residual = probabilities.copy()
    expected_residual.loc[eliminated] = 0.0
    pd.testing.assert_frame_equal(
        result.residual_probabilities, expected_residual
    )

    channel_count = int((metadata.transition_type < 0).sum())
    for cell in range(probabilities.shape[1]):
        original = np.zeros((global_size, global_size))
        channels = np.zeros((global_size, channel_count))
        residual = np.zeros((global_size, global_size + channel_count))
        channel = 0
        for transition, row in metadata.iterrows():
            source = row.source_level_idx
            destination = row.destination_level_idx
            probability = probabilities.iloc[transition, cell]
            retained_probability = result.residual_probabilities.iloc[
                transition, cell
            ]
            if row.transition_type >= 0:
                original[source, destination] += probability
                residual[source, destination] += retained_probability
            else:
                channels[source, channel] += probability
                residual[source, global_size + channel] += retained_probability
                channel += 1
        expected = np.linalg.solve(
            np.eye(len(states)) - original[np.ix_(states, states)],
            channels[states],
        )
        totals = residual.sum(axis=1)
        np.divide(
            residual, totals[:, None], out=residual, where=totals[:, None] != 0
        )
        composed = residual.copy()
        composed[selected] = result.matrix[cell] @ residual[selected]
        actual = np.linalg.solve(
            np.eye(len(states)) - composed[np.ix_(states, states)],
            composed[states, global_size:],
        )
        np.testing.assert_allclose(actual, expected, rtol=2e-14, atol=2e-14)


def test_asymmetric_near_trapping_cycle_preserves_tiny_exit_ratio() -> None:
    """Tiny 1:3 exit rates give asymptotic 1/4:3/4 branching, without cancellation.

    Exact two-state formula retains each actual exit rate; 1e-14 absolute
    tolerance covers final float64 rounding, not an error amplified by trapping.
    """
    first_exit, second_exit = 1e-12, 3e-12
    metadata = pd.DataFrame(
        {
            "source_level_idx": [0, 0, 1, 1],
            "destination_level_idx": [1, -1, 0, -1],
            "transition_type": [1, -1, 1, -1],
            "source": [(1, 0, 0)] * 2 + [(1, 1, 0)] * 2,
        }
    )
    result = create_absorbing_probs(
        pd.DataFrame(
            [1 - first_exit, first_exit, 1 - second_exit, second_exit]
        ),
        metadata,
        mode="partial",
        elements=["H"],
    )
    denominator = first_exit + second_exit - first_exit * second_exit
    expected = (
        np.array(
            [
                [first_exit, (1 - first_exit) * second_exit],
                [(1 - second_exit) * first_exit, second_exit],
            ]
        )
        / denominator
    )
    np.testing.assert_allclose(result.matrix[0], expected, rtol=0, atol=1e-14)
