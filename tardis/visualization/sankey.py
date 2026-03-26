import logging
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple, Union

import pandas as pd
import plotly.graph_objects as go
from collections import defaultdict, Counter


logger = logging.getLogger(__name__)


def create_mock_data() -> pd.DataFrame:
    """Create a mock tracker_full_df-like DataFrame for unit-testing."""
    return pd.DataFrame(
        {
            "packet_id": [0, 0, 0, 0, 1, 1],
            "event_id": [0, 1, 2, 3, 0, 1],
            "interaction_type": [
                "EMISSION",
                "ESCATTERING",
                "BOUNDARY",
                "ESCAPE",
                "EMISSION",
                "LINE",
            ],
            "status": ["IN_PROCESS", "IN_PROCESS", "IN_PROCESS", "EMITTED", "IN_PROCESS", "EMITTED"],
        }
    )


def get_tracking_df(sim_or_df: Union[pd.DataFrame, Any]) -> pd.DataFrame:
    """Normalize input to a tracker_full_df DataFrame."""
    if isinstance(sim_or_df, pd.DataFrame):
        return sim_or_df

    if hasattr(sim_or_df, "transport"):
        transport_state = sim_or_df.transport.transport_state
    elif hasattr(sim_or_df, "transport_state"):
        transport_state = sim_or_df.transport_state
    else:
        raise AttributeError(
            "Simulation object has no transport state or dataframe attribute"
        )

    tracker_full_df = getattr(transport_state, "tracker_full_df", None)
    if tracker_full_df is None:
        raise ValueError(
            "No tracking data found. Enable config.montecarlo.tracking.track_rpacket = True"
        )
    return tracker_full_df.reset_index()


def get_packet_history(tracker_full_df: Union[pd.DataFrame, Any], packet_id: int) -> pd.DataFrame:
    """Return all tracking rows for a single packet by packet_id."""
    df = get_tracking_df(tracker_full_df)
    if "packet_id" not in df.columns:
        raise KeyError("tracker_full_df needs 'packet_id' column")
    packet_df = df[df["packet_id"] == packet_id].sort_values("event_id")
    if packet_df.empty:
        raise ValueError(f"No packet with packet_id={packet_id} found")
    return packet_df


def filter_boundary_interactions(df: pd.DataFrame) -> pd.DataFrame:
    """Drop boundary events from a packet tracking history."""
    return df[df["interaction_type"] != "BOUNDARY"].copy()


def find_packet_with_diverse_interactions(
    tracker_full_df: Union[pd.DataFrame, Any],
    min_unique_interactions: int = 3,
    final_status: Optional[Union[str, Sequence[str]]] = "EMITTED",
    include_interactions: Optional[Sequence[str]] = None,
    exclude_interactions: Optional[Sequence[str]] = None,
    ignore_boundary: bool = True,
) -> int:
    """Select a packet ID with at least `min_unique_interactions` types."""
    df = get_tracking_df(tracker_full_df)

    if ignore_boundary:
        df = filter_boundary_interactions(df)

    if include_interactions is not None:
        df = df[df["interaction_type"].isin(include_interactions)]

    if exclude_interactions is not None:
        df = df[~df["interaction_type"].isin(exclude_interactions)]

    if isinstance(final_status, str):
        final_status_set = {final_status}
    elif final_status is None:
        final_status_set = None
    else:
        final_status_set = set(final_status)

    candidates: List[Tuple[int, int]] = []
    for packet_id, packet in df.groupby("packet_id"):
        history_types = packet["interaction_type"].unique()
        if len(history_types) < min_unique_interactions:
            continue

        if final_status_set is not None:
            last_status = packet.sort_values("event_id")["status"].iloc[-1]
            if last_status not in final_status_set:
                continue

        candidates.append((packet_id, len(history_types)))

    if not candidates:
        raise ValueError(
            "No packet found with the requested diversity/filter constraints"
        )

    # choose packet with maximum unique interactions and fewest events to keep plot manageable
    candidates.sort(key=lambda x: (-x[1], x[0]))
    return candidates[0][0]


def make_sankey_from_packet_history(
    packet_history_df: pd.DataFrame,
    include_interactions: Optional[Sequence[str]] = None,
    exclude_interactions: Optional[Sequence[str]] = None,
    ignore_boundary: bool = True,
    include_start_end: bool = True,
) -> go.Figure:
    """Create a Sankey diagram for a single packet history DataFrame."""
    df = packet_history_df.copy()

    if ignore_boundary:
        df = filter_boundary_interactions(df)

    if include_interactions is not None:
        df = df[df["interaction_type"].isin(include_interactions)]

    if exclude_interactions is not None:
        df = df[~df["interaction_type"].isin(exclude_interactions)]

    if df.empty:
        raise ValueError("Packet history has no interactions after filtering")

    states = [str(s) for s in df["interaction_type"].tolist()]

    # Optionally wrap with START -> ... -> END steps for clarity
    nodes = []
    if include_start_end:
        nodes.append("START")
    for state in states:
        if state not in nodes:
            nodes.append(state)
    if include_start_end:
        nodes.append("END")

    def _node_idx(label: str) -> int:
        return nodes.index(label)

    edge_counts: Dict[Tuple[str, str], int] = {}

    if include_start_end:
        edge_counts[("START", states[0])] = edge_counts.get(("START", states[0]), 0) + 1

    for i in range(len(states) - 1):
        src = states[i]
        tgt = states[i + 1]
        edge_counts[(src, tgt)] = edge_counts.get((src, tgt), 0) + 1

    if include_start_end:
        edge_counts[(states[-1], "END")] = edge_counts.get((states[-1], "END"), 0) + 1

    sources = [ _node_idx(src) for src, _ in edge_counts.keys() ]
    targets = [ _node_idx(tgt) for _, tgt in edge_counts.keys() ]
    values = list(edge_counts.values())

    fig = go.Figure(
        data=[
            go.Sankey(
                node=dict(label=nodes, pad=20, thickness=20, line=dict(color="black", width=0.5)),
                link=dict(source=sources, target=targets, value=values),
            )
        ]
    )

    fig.update_layout(
        title="Packet Interaction History Sankey",
        font=dict(size=12),
    )

    return fig


def plot_packet_history_sankey(
    sim_or_df: Union[pd.DataFrame, Any],
    packet_id: Optional[int] = None,
    min_unique_interactions: int = 3,
    final_status: Optional[Union[str, Sequence[str]]] = "EMITTED",
    include_interactions: Optional[Sequence[str]] = None,
    exclude_interactions: Optional[Sequence[str]] = None,
    ignore_boundary: bool = True,
    include_start_end: bool = True,
) -> go.Figure:
    """
    Generate a Sankey diagram of packet interaction history.

    Works with:
    - simulation object (with tracking enabled)
    - pandas DataFrame (tracker_full_df)

    Must enable: config.montecarlo.tracking.track_rpacket = True

    Parameters
    ----------
    sim_or_df : Simulation or pd.DataFrame
        TARDIS simulation object or tracker_full_df DataFrame.
    packet_id : int, optional
        Specific packet ID to plot. If None, auto-selects diverse packet.
    min_unique_interactions : int, optional
        Minimum unique interactions for auto-selection (default 3).
    final_status : str or list, optional
        Final packet status filter (default "EMITTED").
    include_interactions : list, optional
        Interaction types to include.
    exclude_interactions : list, optional
        Interaction types to exclude.
    ignore_boundary : bool, optional
        Filter out boundary interactions (default True).
    include_start_end : bool, optional
        Add START/END nodes to diagram (default True).

    Returns
    -------
    plotly.graph_objects.Figure
        Sankey diagram figure.

    Examples
    --------
    >>> from tardis.visualization import plot_packet_history_sankey
    >>> fig = plot_packet_history_sankey(sim)
    >>> fig.show()
    """
    """Top-level function to create a packet history Sankey diagram."""

    df = get_tracking_df(sim_or_df)

    if packet_id is None:
        packet_id = find_packet_with_diverse_interactions(
            df,
            min_unique_interactions=min_unique_interactions,
            final_status=final_status,
            include_interactions=include_interactions,
            exclude_interactions=exclude_interactions,
            ignore_boundary=ignore_boundary,
        )

    packet_history_df = get_packet_history(df, packet_id)
    fig = make_sankey_from_packet_history(
        packet_history_df,
        include_interactions=include_interactions,
        exclude_interactions=exclude_interactions,
        ignore_boundary=ignore_boundary,
        include_start_end=include_start_end,
    )
    fig.update_layout(title=f"Packet {packet_id} Interaction History")
    return fig


def print_sequence(df: pd.DataFrame) -> None:
    """Print a linear packet interaction sequence."""
    sequence = df["interaction_type"].astype(str).tolist()
    print(" → ".join(sequence))


def get_line_data_from_simulation(sim_or_df: Union[pd.DataFrame, Any]) -> Optional[pd.DataFrame]:
    """Extract line atomic number and ion number data from simulation."""
    try:
        if hasattr(sim_or_df, "plasma"):
            lines = sim_or_df.plasma.atomic_data.lines
        elif hasattr(sim_or_df, "plasma_solver"):
            lines = sim_or_df.plasma_solver.atomic_data.lines
        else:
            return None
        
        lines_df = lines.reset_index()
        return lines_df[["line_id", "atomic_number", "ion_number"]]
    except Exception:
        return None


def find_carbon_line_ids(lines_df: pd.DataFrame) -> Set[int]:
    """Find all Carbon (atomic number 6) line IDs."""
    if lines_df is None:
        return set()
    
    # Carbon I is atomic_number=6, ion_number=0
    carbon_i = lines_df[
        (lines_df["atomic_number"] == 6) & 
        (lines_df["ion_number"] == 0)
    ]["line_id"].values
    
    return set(carbon_i)


def get_packets_with_carbon_interaction(
    tracker_full_df: pd.DataFrame,
    lines_df: Optional[pd.DataFrame] = None
) -> Dict[int, int]:
    """
    Find packets with Carbon I interactions and return event_id of first C I.
    
    Returns dict: {packet_id: event_id_of_first_carbon_i}
    """
    if lines_df is None:
        return {}
    
    carbon_line_ids = find_carbon_line_ids(lines_df)
    if not carbon_line_ids:
        return {}
    
    df = tracker_full_df.copy()
    if "packet_id" not in df.columns:
        df = df.reset_index()
    
    # Filter for LINE interactions with Carbon I line_emit_id
    df_lines = df[(df["interaction_type"] == "LINE")]
    
    packets_with_c = {}
    for packet_id, packet_group in df_lines.groupby("packet_id"):
        # Find first C I interaction (by line_emit_id)
        c_events = packet_group[packet_group["line_emit_id"].isin(carbon_line_ids)]
        if not c_events.empty:
            first_c_event = c_events.iloc[0]
            packets_with_c[packet_id] = first_c_event["event_id"]
    
    return packets_with_c


def get_post_interaction_flows(
    tracker_full_df: pd.DataFrame,
    packets_with_carbon: Dict[int, int],
    ignore_boundary: bool = True
) -> Dict[Tuple[str, str], int]:
    """
    Build aggregated transition counts from all packets starting after their C I event.
    
    Returns dict: {(from_interaction, to_interaction): packet_count}
    """
    df = tracker_full_df.copy()
    if "packet_id" not in df.columns:
        df = df.reset_index()
    
    transitions = Counter()
    
    for packet_id, carbon_event_id in packets_with_carbon.items():
        # Get packet history after C I
        packet_history = df[df["packet_id"] == packet_id].sort_values("event_id")
        post_carbon = packet_history[packet_history["event_id"] > carbon_event_id]
        
        if post_carbon.empty:
            continue
        
        if ignore_boundary:
            post_carbon = filter_boundary_interactions(post_carbon)
        
        # Build transitions
        interactions = post_carbon["interaction_type"].tolist()
        if len(interactions) > 1:
            for i in range(len(interactions) - 1):
                transitions[(interactions[i], interactions[i + 1])] += 1
        else:
            # Single interaction after C I - path from C I to this interaction
            if interactions:
                transitions[("C I", interactions[0])] += 1
    
    return dict(transitions)


def make_aggregated_sankey_from_flows(
    flows: Dict[Tuple[str, str], int],
    title: str = "Packet Interactions After C I (Carbon Absorption)",
    color_map: Optional[Dict[str, str]] = None
) -> go.Figure:
    """
    Create Sankey diagram from aggregated transition flows.
    
    Parameters
    ----------
    flows : dict
        {(from_state, to_state): count}
    title : str
        Plot title
    color_map : dict, optional
        Custom colors for interaction types
        
    Returns
    -------
    go.Figure
        Sankey diagram
    """
    if not flows:
        raise ValueError("No flows provided")
    
    # Default color palette
    default_colors = {
        "C I": "#1f77b4",
        "LINE": "#ff7f0e",
        "ESCATTERING": "#2ca02c",
        "CONTINUUM_PROCESS": "#d62728",
        "ESCAPE": "#9467bd",
        "EMISSION": "#8c564b",
        "NO_INTERACTION": "#e377c2"
    }
    if color_map:
        default_colors.update(color_map)
    
    # Build node list
    nodes = set()
    for src, tgt in flows.keys():
        nodes.add(src)
        nodes.add(tgt)
    nodes = sorted(list(nodes))
    
    def _node_idx(label: str) -> int:
        return nodes.index(label)
    
    # Build source, target, value lists
    sources = []
    targets = []
    values = []
    
    for (src, tgt), count in flows.items():
        sources.append(_node_idx(src))
        targets.append(_node_idx(tgt))
        values.append(count)
    
    # Get node colors
    node_colors = [default_colors.get(node, "#cccccc") for node in nodes]
    
    # Create Sankey
    fig = go.Figure(
        data=[
            go.Sankey(
                node=dict(
                    label=nodes,
                    color=node_colors,
                    pad=25,
                    thickness=20,
                    line=dict(color="black", width=0.5)
                ),
                link=dict(
                    source=sources,
                    target=targets,
                    value=values,
                    color=[f"rgba(200,200,200,0.4)" for _ in values]
                ),
            )
        ]
    )
    
    fig.update_layout(
        title=title,
        font=dict(size=12),
        height=600,
        width=1000
    )
    
    return fig


def plot_multi_packet_sankey_after_carbon(
    sim_or_df: Union[pd.DataFrame, Any],
    ignore_boundary: bool = True,
) -> go.Figure:
    """
    Create Sankey showing aggregated packet evolution after first Carbon I absorption.
    
    Advanced visualization: shows all packets' paths following their C I interaction,
    with flow thickness proportional to packet count.
    
    Must enable: config.montecarlo.tracking.track_rpacket = True
    
    Parameters
    ----------
    sim_or_df : Simulation or pd.DataFrame
        TARDIS simulation or tracker_full_df DataFrame
    ignore_boundary : bool
        Filter out boundary interactions (default True)
    
    Returns
    -------
    go.Figure
        Multi-packet aggregated Sankey diagram
    
    Examples
    --------
    >>> from tardis.visualization import plot_multi_packet_sankey_after_carbon
    >>> fig = plot_multi_packet_sankey_after_carbon(sim)
    >>> fig.show()
    """
    df = get_tracking_df(sim_or_df)
    
    # Get line data for Carbon identification
    lines_df = None
    if not isinstance(sim_or_df, pd.DataFrame):
        lines_df = get_line_data_from_simulation(sim_or_df)
    
    # Find packets with C I interactions
    packets_with_c = get_packets_with_carbon_interaction(df, lines_df)
    
    if not packets_with_c:
        raise ValueError(
            "No packets with Carbon I interactions found. "
            "Check that tracking is enabled and simulation includes C I interactions."
        )
    
    logger.info(f"Found {len(packets_with_c)} packets with Carbon I interactions")
    
    # Build aggregated flows
    flows = get_post_interaction_flows(df, packets_with_c, ignore_boundary=ignore_boundary)
    
    if not flows:
        raise ValueError("No flows found after Carbon I interactions")
    
    # Create Sankey
    fig = make_aggregated_sankey_from_flows(flows)
    
    return fig



if __name__ == "__main__":
    demo_df = create_mock_data()

    packet_id = find_packet_with_diverse_interactions(demo_df)
    history = get_packet_history(demo_df, packet_id)
    filtered_history = filter_boundary_interactions(history)

    print("Selected packet:", packet_id)
    print_sequence(filtered_history)

    fig = make_sankey_from_packet_history(filtered_history)
    fig.show()   