from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from tardis.transport.montecarlo.packets.trackers.tracker_full_util import (
    generate_tracker_full_list,
    trackers_full_to_df,
    tracker_full_df2tracker_last_interaction_df,
)


@dataclass
class TrackerFullState:
    """
    Holds a list of TrackerFull instances and DataFrames.
    """

    trackers: list
    full_df: pd.DataFrame | None = None
    last_interaction_df: pd.DataFrame | None = None
    has_tracked: bool = False


class TrackerFullSolver:
    """
    Solver for TrackerFull. Holds configuration and generates tracker states.
    """

    def __init__(
        self,
        initial_length: int = 1000,
        extend_factor: int = 2,
    ) -> None:
        self.initial_length = initial_length
        self.extend_factor = extend_factor

    def create_tracker_state(self, num_packets: int) -> TrackerFullState:
        trackers = generate_tracker_full_list(
            num_packets, self.initial_length, self.extend_factor
        )
        return TrackerFullState(trackers)

    def process_tracker_state(
        self, tracker_state: TrackerFullState
    ) -> None:
        """
        Process a completed tracker state and add DataFrames to it.
        Removes the tracker list after processing.
        """
        # Create DataFrames from trackers
        tracker_state.full_df = trackers_full_to_df(tracker_state.trackers)
        tracker_state.last_interaction_df = tracker_full_df2tracker_last_interaction_df(
            tracker_state.full_df
        )

        # Clear trackers and mark as tracked
        tracker_state.trackers.clear()
        tracker_state.has_tracked = True
