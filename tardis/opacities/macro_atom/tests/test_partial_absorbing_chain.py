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
    sources = [0, 0, 2, 2, 2, 4, 4, 4, 5, 5, 6, 6, 6, 9, 12, 12]
    destinations = [2, -1, 4, 9, -1, 0, 5, -1, 6, -1, 0, 12, -1, -1, 6, -1]
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
    # fmt: off
    probabilities = pd.DataFrame({
        10: [0.0, 1.0, 0.3, 0.2, 0.5, 0.4, 0.1, 0.5,
             0.3, 0.7, 0.1, 0.1, 0.8, 1.0, 0.25, 0.75],
        30: [1.0, 0.0, 0.2, 0.1, 0.7, 0.7, 0.1, 0.2,
             0.6, 0.4, 0.05, 0.05, 0.9, 1.0, 0.4, 0.6],
    })
    # fmt: on
    return probabilities, metadata


@pytest.mark.parametrize(
    "mode,elements,include_thermal_pool,include_ionization_pool",
    [
        ("full", [], False, False),
        ("partial", ["H"], False, False),
        ("partial", ["H"], True, False),
        ("partial", ["H"], False, True),
        ("partial", ["H"], True, True),
        ("partial", ["He"], False, False),
        ("partial", ["Fe", "He", "H", "H"], True, True),
        ("explicit", [], False, False),
    ],
)
def test_final_emission_law(
    pool_network: tuple[pd.DataFrame, pd.DataFrame],
    mode: str,
    elements: list[str],
    include_thermal_pool: bool,
    include_ionization_pool: bool,
) -> None:
    """Selection preserves every eventual emission channel in a mixed network.

    Regime: two cells, multiple ions, gapped indices including zero, an
    emission-only destination, zero residual/internal rows, and pool re-entry.
    Either shared pool can be reduced independently; exits from included pools
    to unselected elements and every physical emission remain represented.
    Verify the mask and maps directly, then compare final channel probabilities
    with a dense solve of the original network. The reference never uses the
    reduced matrix. A 2e-14 tolerance allows floating-point LU roundoff.
    """
    probabilities, metadata = pool_network
    result = create_absorbing_probs(
        probabilities,
        metadata,
        mode=mode,
        elements=elements,
        include_thermal_pool=include_thermal_pool,
        include_ionization_pool=include_ionization_pool,
        max_workers=2,
    )
    serial = create_absorbing_probs(
        probabilities,
        metadata,
        mode=mode,
        elements=elements,
        include_thermal_pool=include_thermal_pool,
        include_ionization_pool=include_ionization_pool,
    )
    np.testing.assert_array_equal(result.matrix, serial.matrix)

    states = np.unique(metadata.source_level_idx)
    global_size = int(states.max()) + 1
    atomic_numbers = {"H": 1, "He": 2, "Fe": 26}
    source_identifiers = metadata.source.str[0]
    selected_identifiers: set[int | str] = {
        atomic_numbers[symbol] for symbol in elements
    }
    if include_thermal_pool:
        selected_identifiers.add("k")
    if include_ionization_pool:
        selected_identifiers.add("i")
    if mode == "full":
        selected = np.arange(global_size)
    else:
        selected_rows = source_identifiers.isin(selected_identifiers)
        selected = np.unique(metadata.loc[selected_rows, "source_level_idx"])
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


@pytest.mark.parametrize(
    "mode,include_ionization_pool",
    [
        ("partial", True),
        ("partial", False),
        ("full", False),
        ("explicit", False),
    ],
)
@pytest.mark.parametrize("indirect_rate", [0.6, 0.1])
def test_signed_coefficients_preserve_algebraic_emission_law(
    mode: str,
    include_ionization_pool: bool,
    indirect_rate: float,
) -> None:
    """Preserve full-mode signed-coefficient compatibility in every mode.

    Claim: direct H0->i=-0.1 and indirect H0->H1->i=0.3 combine to
    +0.2; the remaining first-exit masses are 0.5 at H0 and 0.3 at H1.
    Regime: net signed internal coefficients, a shared ionization pool, and
    explicit He re-entry. These are algebraic coefficients, not path weights.
    Verification: exact geometric-series final emissions with a 2e-14
    float64 LU tolerance. Excluding i exposes a negative retained channel;
    reducing the indirect path to 0.05 leaves a negative reduced exit mass.
    Both are retained without rejection or clipping, like full mode. This
    verifies algebraic compatibility, not a categorical sampling law.
    """
    levels = [(1, 0, 0), (1, 0, 1), ("i", -99, -99), (2, 0, 0)]
    sources = [0, 0, 0, 0, 1, 1, 2, 3, 3]
    destinations = [2, 1, 3, -1, 2, -1, -1, 0, -1]
    metadata = pd.DataFrame(
        {
            "source_level_idx": sources,
            "destination_level_idx": destinations,
            "transition_type": [1 if idx >= 0 else -1 for idx in destinations],
            "source": [levels[idx] for idx in sources],
        }
    )
    probabilities = pd.DataFrame(
        [-0.1, indirect_rate, 0.1, 1 - indirect_rate, 0.5, 0.5, 1, 0.25, 0.75]
    )
    options = {
        "mode": mode,
        "elements": ["H"] if mode == "partial" else [],
        "include_ionization_pool": include_ionization_pool,
    }
    result = create_absorbing_probs(probabilities, metadata, **options)
    # Compose the algebraic kernels, then sum explicit He returns. Signed
    # entries are not frequencies that the ordinary CDF sampler can reproduce.
    retained = np.zeros((4, 8))
    emission = 4
    for transition, row in metadata.iterrows():
        destination = row.destination_level_idx
        if destination < 0:
            destination = emission
            emission += 1
        retained[row.source_level_idx, destination] += (
            result.residual_probabilities.iloc[transition, 0]
        )
    retained /= retained.sum(axis=1, keepdims=True)
    selected = result.reduced_to_macro_level_idx
    composed = retained.copy()
    composed[selected] = result.matrix[0] @ retained[selected]
    actual = np.linalg.solve(np.eye(4) - composed[:, :4], composed[:, 4:])
    first_emission = (
        np.array(
            [
                1 - indirect_rate,
                indirect_rate / 2,
                indirect_rate / 2 - 0.1,
                0.075,
            ]
        )
        / 0.975
    )
    expected = np.vstack(
        [
            first_emission,
            [0, 0.5, 0.5, 0],
            [0, 0, 1, 0],
            0.25 * first_emission + [0, 0, 0, 0.75],
        ]
    )
    np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-14)
    if mode == "partial" and include_ionization_pool:
        np.testing.assert_allclose(
            result.matrix[0],
            [
                [
                    1.1 - indirect_rate,
                    indirect_rate / 2,
                    indirect_rate / 2 - 0.1,
                ],
                [0, 0.5, 0.5],
                [0, 0, 1],
            ],
            rtol=0,
            atol=2e-14,
        )
    else:
        # In explicit mode, or partial mode with i outside the selection, the
        # original -0.1 H0->i coefficient must survive unchanged.
        if mode != "full":
            assert result.residual_probabilities.iloc[0, 0] == -0.1
