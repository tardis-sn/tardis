import pandas as pd
import panel as pn
from energy_packet_filtering import filter_by_origin, create_histogram, interactive_plot

# helper function to create a sample df
def create_sample_df():
    data = {
        'nu': [1e14, 2e14, 3e14, 4e14],
        'energy': [10, 20, 30, 40],
        'wavelength': [300e-9, 400e-9, 500e-9, 600e-9],
        'origin_atomic_number': [14, 20, 14, pd.NA]
    }
    return pd.DataFrame(data)


def test_filter_by_origin():
    """
    Test that filter_by_origin returns only rows with the given atomic number.
    """
    df = create_sample_df()
    atomic_number = 14
    filtered_df = filter_by_origin(df, atomic_number)
    # Expected DataFrame is the subset where origin_atomic_number equals 14.
    expected = df[df['origin_atomic_number'] == 14].copy()
    pd.testing.assert_frame_equal(filtered_df, expected)


def test_create_histogram():
    """
    Test that create_histogram returns a Plotly figure with the correct title.
    """
    df = create_sample_df()
    atomic_number = 14
    # Filter the sample DataFrame to mimic expected input.
    df_filtered = filter_by_origin(df, atomic_number)
    fig = create_histogram(df_filtered, atomic_number)
    # Check that the figure title contains the correct atomic number.
    assert f"Energy Packets for Atomic Number {atomic_number}" in fig.layout.title.text

def test_interactive_plot_non_empty(monkeypatch):
    """
    Test interactive_plot when packets exist for the specified atomic number.
    This should return a Plotly pane.
    """
    sample_df = create_sample_df()
    # Patch the global df_energy with sample DataFrame.
    monkeypatch.setattr("energy_packet_filtering.df_energy", sample_df)

    pane = interactive_plot(atomic_number=14)
    # Check that the pane has an 'object' attribute (the Plotly figure).
    assert hasattr(pane, "object")
    fig = pane.object
    # Verify the figure title includes the correct atomic number.
    assert f"Energy Packets for Atomic Number 14" in fig.layout.title.text
