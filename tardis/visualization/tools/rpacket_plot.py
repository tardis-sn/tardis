import logging
from typing import Dict, List, Optional, Tuple

import astropy.units as u
import numpy as np
import pandas as pd
import plotly.graph_objects as go

from tardis.transport.montecarlo.packets.radiative_packet import InteractionType


class RPacketPlotter:
    """
    Animated plotting interface for Monte Carlo packet trajectories.

    This class creates an animated plot using Plotly to visualize the trajectories
    of real packets as they travel through the ejecta starting from the photosphere.
    The visualization shows packet interactions and shell structures in 2D space.

    Parameters
    ----------
    sim : Simulation
        TARDIS simulation object containing transport state and configuration.
    no_of_packets : int
        Number of packets to include in the visualization.

    Attributes
    ----------
    no_of_packets : int
        Number of packets being visualized.
    sim : Simulation
        Reference to the TARDIS simulation object.
    interaction_from_num : dict[int, dict[str, str | float]]
        Mapping from interaction type integers to display properties.
    theme_colors : dict[str, dict[str, str]]
        Color schemes for light and dark themes.
    play_button : dict
        Plotly button configuration for animation play control.
    pause_button : dict
        Plotly button configuration for animation pause control.
    fig : go.Figure, optional
        The generated Plotly figure object.
    """

    def __init__(self, sim, no_of_packets: int) -> None:
        """
        Initialize the RPacket Plotter.

        Parameters
        ----------
        sim : Simulation
            TARDIS simulation object generated using the run_tardis function.
        no_of_packets : int
            Number of packets to be used for plotting. Must be positive.

        Raises
        ------
        ValueError
            If no_of_packets is not positive.
        """
        if no_of_packets <= 0:
            msg = "no_of_packets must be positive"
            raise ValueError(msg)

        self.no_of_packets = no_of_packets
        self.sim = sim
        # Create a dictionary mapping interaction type string values to display properties
        # This allows direct lookup using string interaction_type values
        self.interaction_from_num = {
            "NO_INTERACTION": {
                "text": "No Interaction",
                "color": "#2E86AB",
                "color_dark": "#2E86AB",
                "opacity": 0.8,
            },
            "BOUNDARY": {
                "text": "Boundary",
                "color": "#A23B72",
                "color_dark": "#A23B72",
                "opacity": 0.8,
            },
            "LINE": {
                "text": "Line Interaction",
                "color": "#F18F01",
                "color_dark": "#F18F01",
                "opacity": 0.8,
            },
            "ESCATTERING": {
                "text": "E-Scattering",
                "color": "#C73E1D",
                "color_dark": "#C73E1D",
                "opacity": 0.8,
            },
            "CONTINUUM_PROCESS": {
                "text": "Continuum",
                "color": "#6A4C93",
                "color_dark": "#6A4C93",
                "opacity": 0.8,
            },

        }
        self.theme_colors = dict(
            light=dict(
                linecolor="#555",
                gridcolor="#fafafa",
                zerolinecolor="#fafafa",
                color="#000",
                photosphere_line_color="black",
                photosphere_fillcolor="darkgrey",
                shells_line_color="black",
                packet_line_color="darkslategrey",
                plot_bgcolor="#fafafa",
                paper_bgcolor="#fafafa",
                title_font_color="#444",
                legendgrouptitle_color="#444",
                button_bgcolor="#fafafa",
                button_font_color="#2A3F5F",
                slider_font_color="#2A3F5F",
                bordercolor="#BEC8D9",
                slider_bgcolor="#F8FAFC",
                slider_activebgcolor="#DBDDE0",
                slider_currentvalue_color="#2A3F5F",
                font_color="#000",
            ),
            dark=dict(
                linecolor="#050505",
                gridcolor="#111",
                zerolinecolor="#111",
                color="#fafafa",
                photosphere_line_color="#222",
                photosphere_fillcolor="#222",
                shells_line_color="#555",
                packet_line_color="#888",
                plot_bgcolor="#000",
                paper_bgcolor="#000",
                title_font_color="#ccc",
                legendgrouptitle_color="#ccc",
                button_bgcolor="#282828",
                button_font_color="#888",
                slider_font_color="#888",
                bordercolor="black",
                slider_bgcolor="#888",
                slider_activebgcolor="#fafafa",
                slider_currentvalue_color="#fff",
                font_color="#fafafa",
            ),
        )

        self.play_button = {
            "args": [
                None,
                {
                    "frame": {"duration": 500, "redraw": False},
                    "fromcurrent": True,
                    "transition": {
                        "duration": 300,
                        "easing": "quadratic-in-out",
                    },
                },
            ],
            "label": "Play",
            "method": "animate",
        }

        self.pause_button = {
            "args": [
                [None],
                {
                    "frame": {"duration": 0, "redraw": False},
                    "mode": "immediate",
                    "transition": {"duration": 0},
                },
            ],
            "label": "Pause",
            "method": "animate",
        }

    @classmethod
    def from_simulation(
        cls,
        sim,
        no_of_packets: int = 15,
    ) -> "RPacketPlotter":
        """
        Create an RPacketPlotter instance from a TARDIS simulation.

        This factory method creates a plotter instance with validation
        to ensure the simulation has the necessary tracker data.

        Parameters
        ----------
        sim : Simulation
            TARDIS simulation object generated using run_tardis function.
            Must have rpacket tracking enabled.
        no_of_packets : int, optional
            Number of packets to visualize, by default 15.
            If greater than available packets, all available packets are used.

        Returns
        -------
        RPacketPlotter
            Configured plotter instance ready for visualization.

        Raises
        ------
        AttributeError
            If the simulation does not have rpacket tracking enabled.

        Warns
        -----
        UserWarning
            If requested no_of_packets exceeds available packets.
        """
        logger = logging.getLogger(__name__)
        if hasattr(sim.transport.transport_state, "tracker_full_df"):
            if sim.last_no_of_packets >= no_of_packets:
                return cls(sim, no_of_packets)
            logger.warning(
                "no_of_packets specified are more than the actual no of packets "
                "in the model. Using all packets in the model."
            )
            return cls(sim, sim.last_no_of_packets)
        raise AttributeError(
            "There is no attribute named rpacket_tracker in the simulation object "
            "passed. Try enabling the rpacket tracking in the configuration. "
            "To enable rpacket tracking see: "
            "https://tardis-sn.github.io/tardis/io/output/rpacket_tracking.html"
            "#How-to-Setup-the-Tracking-for-the-RPackets?"
        )

    def generate_plot(self, theme: str = "light") -> go.Figure:
        """
        Create an animated plotly plot showing Monte Carlo packet trajectories.

        This method generates a comprehensive visualization showing packet paths,
        interaction types, shell boundaries, and photosphere. The plot includes
        animation controls and customizable themes.

        Parameters
        ----------
        theme : str, optional
            Visual theme for the plot, by default "light".
            Must be either "light" or "dark".

        Returns
        -------
        go.Figure
            Plotly figure object containing the animated packet trajectory plot
            with shells, photosphere, and interaction visualization.

        Raises
        ------
        ValueError
            If theme is not "light" or "dark".
        """
        if theme not in ("light", "dark"):
            msg = f"Theme must be 'light' or 'dark', got '{theme}'"
            raise ValueError(msg)
        self.fig = go.Figure()

        # getting velocity of different shells
        v_shells = self.sim.simulation_state.velocity.to_value(u.km / u.s)

        # getting coordinates and interactions of all packets
        (
            rpacket_x,
            rpacket_y,
            rpacket_interactions,
        ) = self.get_coordinates_multiple_packets(
            self.sim.transport.transport_state.tracker_full_df.loc[
                0 : (self.no_of_packets)
            ],
        )

        # making the coordinate arrays of all packets equal
        (
            rpacket_x,
            rpacket_y,
            rpacket_interactions,
            rpacket_array_max_size,
        ) = self.get_equal_array_size(
            rpacket_x, rpacket_y, rpacket_interactions
        )

        axis_props = dict(
            range=[-1.1 * v_shells[-1], 1.1 * v_shells[-1]],
            title="Velocity (km/s)",
            exponentformat="none",
            color=self.theme_colors[theme]["color"],
            linecolor=self.theme_colors[theme]["linecolor"],
            gridcolor=self.theme_colors[theme]["gridcolor"],
            zerolinecolor=self.theme_colors[theme]["zerolinecolor"],
        )
        # Set axes properties
        self.fig.update_xaxes(scaleanchor="y", scaleratio=1, **axis_props)
        self.fig.update_yaxes(**axis_props)

        # adding the shells and photosphere
        for shell_no in range(len(self.sim.simulation_state.radius.value)):
            shape_props = dict(
                type="circle",
                xref="x",
                yref="y",
                x0=-1 * v_shells[shell_no],
                y0=-1 * v_shells[shell_no],
                x1=v_shells[shell_no],
                y1=v_shells[shell_no],
            )
            if shell_no == 0:
                # photosphere
                self.fig.add_shape(
                    **shape_props,
                    line_color=self.theme_colors[theme][
                        "photosphere_line_color"
                    ],
                    fillcolor=self.theme_colors[theme]["photosphere_fillcolor"],
                    opacity=1,
                )
            elif shell_no == (len(self.sim.simulation_state.radius.value) - 1):
                # outermost shell
                self.fig.add_shape(
                    **shape_props,
                    line_color=self.theme_colors[theme]["shells_line_color"],
                    opacity=1,
                )
            else:
                # remaining shells
                self.fig.add_shape(
                    **shape_props,
                    line_color=self.theme_colors[theme]["shells_line_color"],
                    opacity=0.2,
                )

        # Adding packet trajectory
        for packet_no in range(len(rpacket_x)):
            self.fig.add_trace(
                self.create_packet_scatter(
                    packet_no, rpacket_x, rpacket_y, rpacket_interactions, theme
                )
            )

        # adding legends
        legend_props = dict(
            x=[9999999],
            y=[0],
            legendgroup="a",
            opacity=1,
            mode="lines+markers",
        )
        self.fig.add_trace(
            go.Scatter(
                **legend_props,
                legendgrouptitle=dict(
                    font=dict(
                        color=self.theme_colors[theme]["legendgrouptitle_color"]
                    ),
                    text="Interaction Type:",
                ),
                name="e-scattering",
                hoverlabel=dict(font=dict(color="#222")),
                marker=dict(color=self.interaction_from_num["ESCATTERING"]["color"]),
            )
        )
        self.fig.add_trace(
            go.Scatter(
                **legend_props,
                name="Line Interaction",
                marker=dict(color=self.interaction_from_num["LINE"]["color"]),
            )
        )

        # Set figure size
        self.fig.layout.plot_bgcolor = self.theme_colors[theme]["plot_bgcolor"]
        self.fig.layout.paper_bgcolor = self.theme_colors[theme][
            "paper_bgcolor"
        ]

        self.fig.update_layout(
            width=820,
            height=680,
            title="Packet Trajectories",
            title_font_color=self.theme_colors[theme]["title_font_color"],
            font_color=self.theme_colors[theme]["font_color"],
            updatemenus=[
                dict(
                    type="buttons",
                    xanchor="right",
                    x=0.1,
                    y=0,
                    yanchor="top",
                    direction="left",
                    pad={"r": 10, "t": 87},
                    showactive=False,
                    bgcolor=self.theme_colors[theme]["button_bgcolor"],
                    bordercolor=self.theme_colors[theme]["bordercolor"],
                    font={
                        "color": self.theme_colors[theme]["button_font_color"]
                    },
                    buttons=[self.play_button, self.pause_button],
                )
            ],
        )

        # adding frames
        self.fig.frames = [
            go.Frame(
                data=self.get_frames(
                    frame, rpacket_x, rpacket_y, rpacket_interactions, theme
                ),
                name=frame,
            )
            for frame in range(rpacket_array_max_size + 1)
        ]

        # adding timeline slider
        self.fig.layout.sliders = [
            {
                "active": 0,
                "activebgcolor": self.theme_colors[theme][
                    "slider_activebgcolor"
                ],
                "bgcolor": self.theme_colors[theme]["slider_bgcolor"],
                "bordercolor": self.theme_colors[theme]["bordercolor"],
                "yanchor": "top",
                "xanchor": "left",
                "currentvalue": {
                    "font": {
                        "color": self.theme_colors[theme][
                            "slider_currentvalue_color"
                        ],
                    },
                    "prefix": "Step:",
                    "visible": True,
                    "xanchor": "right",
                },
                "font": {
                    "color": self.theme_colors[theme]["slider_font_color"]
                },
                "transition": {"duration": 300, "easing": "cubic-in-out"},
                "pad": {"b": 10, "t": 50},
                "len": 0.9,
                "x": 0.1,
                "y": 0,
                "steps": self.get_slider_steps(rpacket_array_max_size),
            }
        ]

        return self.fig

    def get_coordinates_with_theta_init(
        self,
        r_track: pd.Series,
        mu_track: pd.Series,
        time: float,
        last_interaction_type: pd.Series,
        theta_initial: float = 0,
    ) -> Tuple[np.ndarray, np.ndarray, List[int]]:
        """
        Generate 2D coordinates for a single packet trajectory.

        Calculates x, y coordinates from radial position and direction cosine
        data, converting from spherical to Cartesian coordinates for visualization.

        Parameters
        ----------
        r_track : pd.Series
            Radial position of packet at each step (in cm).
        mu_track : pd.Series
            Cosine of radial angle of the packet at each step.
        time : float
            Time since explosion occurrence (in seconds).
        last_interaction_type : pd.Series
            Interaction type of the packet at each step.
        theta_initial : float, optional
            Initial launch angle from x-axis at photosphere, by default 0.

        Returns
        -------
        tuple[np.ndarray, np.ndarray, list[int]]
            Three-element tuple containing:
            - x coordinates of the packet at different steps (km/s)
            - y coordinates of the packet at different steps (km/s)
            - interaction types occurring at different points

        Notes
        -----
        Coordinates are converted to velocity units (km/s) by dividing by time
        and applying unit conversion. The trajectory calculation follows the
        spherical geometry described in the TARDIS documentation.
        """
        theta, rpacket_interactions = [], []

        # getting thetas at different steps of the packet movement

        for step_no in range(0, len(r_track)):
            # for the first step the packet is at photosphere, so theta will be equal to the intial angle we are launching the packet from
            if step_no == 0:
                theta.append(theta_initial)
            # for further steps we calculate thetas with the formula derived in the documentation
            # https://tardis-sn.github.io/tardis/analyzing_tardis/visualization/tutorial_montecarlo_packet_visualization.html#Getting-packet-coordinates
            else:
                curr_r = r_track[step_no]
                prev_r = r_track[step_no - 1]
                prev_mu = mu_track[step_no - 1]

                acos_mu = np.acos(prev_mu)
                sin_term = prev_r * np.sin(acos_mu) / curr_r
                new_theta = theta[-1] + acos_mu

                if curr_r < prev_r:
                    new_theta = new_theta - np.pi + np.asin(sin_term)
                else:
                    new_theta += np.asin(-1 * sin_term)
                theta.append(new_theta)


        # converting the thetas into x and y coordinates using radius as radius*cos(theta) and radius*sin(theta) respectively
        rpacket_x = (np.array(r_track)) * np.cos(np.array(theta)) * 1e-5 / time
        rpacket_y = (np.array(r_track)) * np.sin(np.array(theta)) * 1e-5 / time

        # adding interactions at different steps
        # using the change of slope of the trajectory line at different steps, we determine if an interactions happened or not.

        for step_no in range(len(r_track)):
            # when packet is at its starting and ending point in its trajectory, we consider it as no interaction
            if step_no == 0 or step_no == len(r_track) - 1:
                rpacket_interactions.append("NO_INTERACTION")
            else:
                # current slope is the slope of line from previous position of the packet to the current position
                rpacket_interactions.append(last_interaction_type[step_no])

        return rpacket_x, rpacket_y, rpacket_interactions

    def get_coordinates_multiple_packets(
        self, r_packet_tracker: pd.DataFrame
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Generate coordinates for multiple packet trajectories.

        Creates uniformly distributed launch angles and calculates coordinates
        for multiple packets to visualize their collective behavior.

        Parameters
        ----------
        r_packet_tracker : pd.DataFrame
            DataFrame containing packet tracking data with columns for radius,
            direction cosine, and interaction types.

        Returns
        -------
        tuple[np.ndarray, np.ndarray, np.ndarray]
            Three-element tuple containing:
            - Array of x coordinate arrays for multiple packets (object dtype)
            - Array of y coordinate arrays for multiple packets (object dtype)
            - Array of interaction type arrays for multiple packets (object dtype)

        Notes
        -----
        Launch angles are distributed uniformly around the photosphere to
        provide representative coverage of packet trajectories.
        """
        # Remove rows with any NaN values
        r_packet_tracker = r_packet_tracker.dropna()
        
        # for plotting packets at equal intervals throught the circle, we choose thetas distributed uniformly
        thetas = np.linspace(0, 2 * np.pi, self.no_of_packets + 1)
        all_rpackets_x_coords = []
        all_rpackets_y_coords = []
        all_rpackets_interactions_coords = []

        # getting coordinates and interaction arrays for all packets
        for packet_no in range(self.no_of_packets):
            packet_data = r_packet_tracker.loc[packet_no]
            interaction_types = packet_data["interaction_type"]
            mu_data = packet_data["after_mu"]
            
            (
                rpacket_x,
                rpacket_y,
                rpacket_interactions,
            ) = self.get_coordinates_with_theta_init(
                packet_data["radius"],
                mu_data,
                self.sim.simulation_state.time_explosion.value,
                interaction_types,
                thetas[packet_no],
            )
            all_rpackets_x_coords.append(rpacket_x)
            all_rpackets_y_coords.append(rpacket_y)
            all_rpackets_interactions_coords.append(rpacket_interactions)
        return (
            np.array(all_rpackets_x_coords, dtype="object"),
            np.array(all_rpackets_y_coords, dtype="object"),
            np.array(all_rpackets_interactions_coords, dtype="object"),
        )

    def get_equal_array_size(
        self,
        rpacket_x: np.ndarray,
        rpacket_y: np.ndarray,
        interactions: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
        """
        Normalize coordinate arrays to equal size for animation frames.

        Pads shorter trajectories with their final values to match the longest
        trajectory length, enabling synchronized animation across all packets.

        Parameters
        ----------
        rpacket_x : np.ndarray
            Array of x coordinate arrays for different packets.
        rpacket_y : np.ndarray
            Array of y coordinate arrays for different packets.
        interactions : np.ndarray
            Array of interaction type arrays for different packets.

        Returns
        -------
        tuple[np.ndarray, np.ndarray, np.ndarray, int]
            Four-element tuple containing:
            - Normalized x coordinate array (all sub-arrays same length)
            - Normalized y coordinate array (all sub-arrays same length)
            - Normalized interaction types array (all sub-arrays same length)
            - Maximum array size among all packets

        Notes
        -----
        Padding uses the final coordinate values to maintain trajectory continuity
        in the animation visualization.
        """
        rpacket_step_no_array_max_size = max(list(map(len, rpacket_x)))

        for packet_no in range(len(rpacket_x)):
            padding_size = rpacket_step_no_array_max_size - len(
                rpacket_x[packet_no]
            )
            if padding_size > 0:
                x_padding = np.ones(padding_size) * rpacket_x[packet_no][-1]
                y_padding = np.ones(padding_size) * rpacket_y[packet_no][-1]
                interactions_padding = np.full(padding_size, interactions[packet_no][-1])

                rpacket_x[packet_no] = np.append(
                    rpacket_x[packet_no], x_padding
                )
                rpacket_y[packet_no] = np.append(
                    rpacket_y[packet_no], y_padding
                )
                interactions[packet_no] = np.append(
                    interactions[packet_no], interactions_padding
                )
        return (
            rpacket_x,
            rpacket_y,
            interactions,
            rpacket_step_no_array_max_size,
        )

    def get_frames(
        self,
        frame: int,
        rpacket_x: np.ndarray,
        rpacket_y: np.ndarray,
        interactions: np.ndarray,
        theme: str,
    ) -> List[go.Scatter]:
        """
        Create individual animation frames with scatter plot objects.

        Generates a list of scatter plot objects for each packet trajectory
        up to the specified frame number for animation visualization.

        Parameters
        ----------
        frame : int
            Current frame number for trajectory slicing.
        rpacket_x : np.ndarray
            Array of x coordinate arrays for different packets.
        rpacket_y : np.ndarray
            Array of y coordinate arrays for different packets.
        interactions : np.ndarray
            Array of interaction type arrays for different packets.
        theme : str
            Theme name for plot styling.

        Returns
        -------
        list[go.Scatter]
            List of Plotly scatter objects for the current frame.

        Notes
        -----
        Each scatter object represents one packet's trajectory up to the
        current frame, enabling smooth animation transitions.
        """
        return [
            self.create_packet_scatter(
                packet_no, rpacket_x, rpacket_y, interactions, theme, frame=frame
            )
            for packet_no in range(len(rpacket_x))
        ]

    def create_packet_scatter(
        self,
        packet_no: int,
        rpacket_x: np.ndarray,
        rpacket_y: np.ndarray,
        interactions: np.ndarray,
        theme: str,
        frame: Optional[int] = None,
    ) -> go.Scatter:
        """
        Create a scatter plot object for a single packet trajectory.

        Generates a Plotly scatter plot with interaction-based styling for
        visualizing packet trajectories with optional frame-based slicing.

        Parameters
        ----------
        packet_no : int
            Index of the packet to create scatter plot for.
        rpacket_x : np.ndarray
            Array of x coordinate arrays for different packets.
        rpacket_y : np.ndarray
            Array of y coordinate arrays for different packets.
        interactions : np.ndarray
            Array of interaction type arrays for different packets.
        theme : str
            Theme name for plot styling.
        frame : int, optional
            Frame number for trajectory slicing. If None, shows full trajectory.

        Returns
        -------
        go.Scatter
            Plotly scatter object with styled trajectory and interaction markers.

        Notes
        -----
        Marker colors and opacity are determined by interaction types.
        Hover information includes coordinates and last interaction type.
        """
        x_vals = rpacket_x[packet_no]
        y_vals = rpacket_y[packet_no]
        interaction_vals = interactions[packet_no]

        # If frame is specified, slice data
        if frame is not None:
            x_vals = x_vals[:frame]
            y_vals = y_vals[:frame]
            interaction_vals = interaction_vals[:frame]

        text_labels = [
            self.interaction_from_num[interaction]["text"]
            for interaction in interaction_vals
        ]
        opacity_values = [
            self.interaction_from_num[interaction]["opacity"]
            for interaction in interaction_vals
        ]
        color_values = [
            self.interaction_from_num[interaction]["color"]
            for interaction in interaction_vals
        ]

        return go.Scatter(
            x=x_vals,
            y=y_vals,
            mode="markers+lines",
            name=f"Packet {packet_no + 1}",
            showlegend=False,
            hovertemplate=(
                "<b>X</b>: %{x}"
                "<br><b>Y</b>: %{y}<br>"
                "<b>Last Interaction: %{text}</b>"
            ),
            text=text_labels,
            line=dict(color=self.theme_colors[theme]["packet_line_color"]),
            marker=dict(
                opacity=opacity_values,
                color=color_values,
            ),
        )

    def get_slider_steps(self, rpacket_max_array_size: int) -> List[Dict]:
        """
        Generate timeline slider steps for animated plot frames.

        Creates slider step configurations for frame-by-frame animation
        control with consistent transition timing.

        Parameters
        ----------
        rpacket_max_array_size : int
            Maximum size of coordinate array among all packets.

        Returns
        -------
        list[dict]
            List of slider step dictionaries with frame animation arguments.

        Notes
        -----
        Each step includes frame duration, transition settings, and label
        for user interaction with the animation timeline.
        """
        slider_steps = []
        base_step_args = {
            "frame": {"duration": 300, "redraw": False},
            "mode": "immediate",
            "transition": {"duration": 300},
        }
        for step_no in range(rpacket_max_array_size):
            slider_steps.append(
                {
                    "args": [[step_no], base_step_args],
                    "label": step_no,
                    "method": "animate",
                }
            )

        return slider_steps
