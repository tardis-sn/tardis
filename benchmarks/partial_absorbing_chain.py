"""Remote ASV setup/solve and retained-memory comparison on a fixed network."""

import numpy as np
import pandas as pd

from tardis.opacities.macro_atom.absorbing_markov_chain import (
    create_absorbing_probs,
)


class BenchmarkPartialAbsorbingChain:
    """Keep physics and cells fixed while changing the reduction boundary."""

    params = ["full", "H", "He", "H+He", "explicit"]
    param_names = ["selection"]

    def setup(self, selection: str) -> None:
        """Build H/He/Fe levels plus shared i/k pools in eight cells."""
        states = [
            (atomic_number, 0, level)
            for atomic_number in (1, 2, 26)
            for level in range(30)
        ]
        states += [("i", -99, -99), ("k", -99, -99)]
        sources = np.repeat(np.arange(len(states)), 3)
        self.metadata = pd.DataFrame(
            {
                "source": [states[source] for source in sources],
                "source_level_idx": sources,
                "destination_level_idx": np.column_stack(
                    (
                        (np.arange(len(states)) + 1) % len(states),
                        np.full(len(states), 90),
                        np.full(len(states), -1),
                    )
                ).ravel(),
                "transition_type": np.tile([1, 2, -1], len(states)),
            }
        )
        self.probabilities = pd.DataFrame(
            np.tile([0.4, 0.2, 0.4], len(states))[:, None] * np.ones((1, 8))
        )
        self.mode = (
            selection if selection in ("full", "explicit") else "partial"
        )
        self.elements = (
            () if self.mode != "partial" else tuple(selection.split("+"))
        )
        self.result = create_absorbing_probs(
            self.probabilities,
            self.metadata,
            mode=self.mode,
            elements=self.elements,
        )

    def time_create_absorbing_probs(self, selection: str) -> None:
        """Measure full mapping construction, diagnostics, and per-cell LU solve."""
        create_absorbing_probs(
            self.probabilities,
            self.metadata,
            mode=self.mode,
            elements=self.elements,
        )

    def track_selected_states(self, selection: str) -> int:
        """Report selected real/pool state count for interpreting solve cost."""
        return len(self.result.reduced_to_macro_level_idx)

    def track_dense_matrix_bytes(self, selection: str) -> int:
        """Report retained matrix memory summed over the same eight cells."""
        return self.result.matrix.nbytes
