import pandas as pd
import plotly.express as px
import panel as pn
from tardis import run_tardis
from tardis.io.configuration.config_reader import Configuration
from tardis.io.atom_data import download_atom_data
import astropy.units as u

pn.extension("plotly")

def setup_simulation(yaml_file: str = "tardis_example.yml"):
    """runs the TARDIS simulation and prepares dataframes."""
    download_atom_data('kurucz_cd23_chianti_H_He')
    config = Configuration.from_yaml(yaml_file)
    config["montecarlo"]["tracking"]["track_rpacket"] = True
    sim = run_tardis(config, show_progress_bars=False)
    tracker_df = sim.transport.transport_state.rpacket_tracker_df
    df_energy = tracker_df[['nu', 'energy']].copy()
    df_energy['wavelength'] = (df_energy['nu'].values * u.Hz).to(u.m, equivalencies=u.spectral()).value
    df_energy['origin_atomic_number'] = pd.NA

    lines_df = sim.plasma.atomic_data.lines.reset_index()
    interaction_ids = sim.transport.transport_state.last_line_interaction_in_id

    for idx, line_id in enumerate(interaction_ids):
        if line_id != -1:
            try:
                df_energy.at[idx, 'origin_atomic_number'] = lines_df.iloc[line_id]['atomic_number']
            except IndexError:
                df_energy.at[idx, 'origin_atomic_number'] = pd.NA

    return df_energy


def filter_by_origin(df, atomic_number):
    """Filters the energy packet df for a specific atomic number."""
    return df[df['origin_atomic_number'] == atomic_number].copy()


def create_histogram(df_filtered, atomic_number):
    """Creates an interactive Plotly histogram."""
    df_filtered['wavelength_A'] = df_filtered['wavelength'] * 1e10
    print(df_filtered['wavelength_A'])
    fig = px.histogram(
        df_filtered,
        x='wavelength_A',
        nbins=100,
        title=f"Energy Packets for Atomic Number {atomic_number}",
        labels={'wavelength_A': 'Wavelength [Å]'},
        color_discrete_sequence=['crimson']
    )
    fig.update_layout(
        xaxis_title='Wavelength [Å]',
        yaxis_title='Number of Packets',
        bargap=0.05,
    )
    return fig


# Global simulation data for interactive use
df_energy = setup_simulation()

def interactive_plot(atomic_number):
    """Returns a Panel pane: histogram or no data message."""
    filtered_data = filter_by_origin(df_energy, atomic_number)
    if filtered_data.empty:
        return pn.pane.Markdown(f"### No packets found for atomic number {atomic_number}")
    return pn.pane.Plotly(create_histogram(filtered_data, atomic_number), config={'responsive': True})


def build_dashboard():
    """Creates the interactive Panel dashboard."""
    unique_atomic_numbers = sorted(df_energy['origin_atomic_number'].dropna().unique())
    selector = pn.widgets.Select(name="Atomic Number", options=unique_atomic_numbers, value=14)

    dashboard = pn.Column(
        "# Energy Packet Histogram by Wavelength",
        "Select an atomic number to visualize:",
        selector,
        pn.bind(interactive_plot, atomic_number=selector),
    )
    return dashboard

# Only expose dashboard when running as script/server
if __name__.startswith("bokeh"):
    build_dashboard().servable()
