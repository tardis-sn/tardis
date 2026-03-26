import pandas as pd

from tardis.visualization.sankey import (
    create_mock_data,
    filter_boundary_interactions,
    find_packet_with_diverse_interactions,
    get_packet_history,
    make_sankey_from_packet_history,
    plot_packet_history_sankey,
    get_post_interaction_flows,
    make_aggregated_sankey_from_flows,
    find_carbon_line_ids,
)


def test_get_packet_history_and_filter_boundary():
    df = create_mock_data()
    packet_df = get_packet_history(df, 0)
    assert len(packet_df) == 4

    filtered = filter_boundary_interactions(packet_df)
    assert "BOUNDARY" not in filtered["interaction_type"].values
    assert list(filtered["interaction_type"]) == ["EMISSION", "ESCATTERING", "ESCAPE"]


def test_find_packet_with_diverse_interactions():
    df = create_mock_data()
    packet_id = find_packet_with_diverse_interactions(df, min_unique_interactions=3)
    assert packet_id == 0


def test_make_sankey_from_packet_history():
    df = create_mock_data()
    packet_df = get_packet_history(df, 0)
    fig = make_sankey_from_packet_history(packet_df)

    assert fig is not None
    assert fig.data
    sankey = fig.data[0]
    assert sankey.type == "sankey"
    assert len(sankey.node.label) >= 4


def test_plot_packet_history_sankey_auto_find():
    df = create_mock_data()
    fig = plot_packet_history_sankey(df, packet_id=None, min_unique_interactions=3)
    assert fig is not None
    assert fig.data


def test_get_post_interaction_flows():
    """Test multi-packet flow aggregation."""
    # Create test data with multiple packets
    df = pd.DataFrame({
        "packet_id": [0, 0, 0, 0, 1, 1, 1, 2, 2, 2],
        "event_id": [0, 1, 2, 3, 0, 1, 2, 0, 1, 2],
        "interaction_type": [
            "C I", "LINE", "ESCATTERING", "ESCAPE",
            "C I", "ESCATTERING", "ESCAPE",
            "C I", "ESCATTERING", "ESCAPE"
        ],
        "line_emit_id": [1, -1, -1, -1, 1, -1, -1, 1, -1, -1]
    })
    
    # Simulate packets with C I at event 0
    packets_with_c = {0: 0, 1: 0, 2: 0}
    
    flows = get_post_interaction_flows(df, packets_with_c, ignore_boundary=True)
    
    # Should have transitions: C I -> LINE, C I -> ESCATTERING, etc.
    assert flows is not None
    assert len(flows) > 0


def test_make_aggregated_sankey_from_flows():
    """Test Sankey creation from aggregated flows."""
    flows = {
        ("C I", "LINE"): 5,
        ("C I", "ESCATTERING"): 3,
        ("LINE", "ESCAPE"): 5,
        ("ESCATTERING", "ESCAPE"): 3,
    }
    
    fig = make_aggregated_sankey_from_flows(flows)
    
    assert fig is not None
    assert fig.data
    sankey = fig.data[0]
    assert sankey.type == "sankey"
    # Should have at least 4 nodes (C I, LINE, ESCATTERING, ESCAPE)
    assert len(sankey.node.label) >= 4
    # Should have 4 links
    assert len(sankey.link.source) == 4


def test_find_carbon_line_ids():
    """Test Carbon line identification."""
    lines_df = pd.DataFrame({
        "line_id": [1, 2, 3, 4],
        "atomic_number": [6, 6, 8, 6],
        "ion_number": [0, 1, 0, 2]
    })
    
    c_ids = find_carbon_line_ids(lines_df)
    assert 1 in c_ids  # Carbon I
    assert 2 not in c_ids  # Carbon II
    assert 3 not in c_ids  # Oxygen I

