import logging
import math

import astropy.units as u
import numpy as np
import plotly.graph_objects as go


class RPacketPlotter:
    """
    Plotting interface for the plotting Montecarlo packets. It creates an animated plot using plotly for the
    trajectories of the real packets as they travel through the ejecta starting from the photosphere.
    """

    def __init__(self, sim, no_of_packets):
        """
        Initializes the RPacket Plotter using the simulation object generated using the run_tardis function.

        Parameters
        ----------
        sim : tardis.simulation.Simulation
            simulation object generated using the run_tardis function.
        no_of_packets : int
            number of packets to be used for plotting.
        """
        self.no_of_packets = no_of_packets
        self.sim = sim
        self.interaction_from_num = [
            {"text": "No Interaction", "color": "darkslategrey", "opacity": 0},
            {"text": "e-Scattering", "color": "#3366FF", "opacity": 1},
            {"text": "Line Interaction", "color": "#FF3300", "opacity": 1},
        ]
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
    def from_simulation(cls, sim, no_of_packets=15):
        """
        Creates an instance of RPacketPlotter from a TARDIS simulation object.

        Parameters
        ----------
        sim : tardis.simulation.Simulation
            TARDIS Simulation object generated using run_tardis function.
        no_of_packets : int
            number of packets to be used for plotting.

        Returns
        -------
        RPacketPlotter
        """
        logger = logging.getLogger(__name__)
        if hasattr(sim.transport.transport_state, "rpacket_tracker_df"):
            if sim.last_no_of_packets >= no_of_packets:
                return cls(sim, no_of_packets)
            logger.warning(
                """
                no_of_packets specified are more than the actual no of packets in the model. Using all packets in the model.
                """
            )
            return cls(sim, sim.last_no_of_packets)
        raise AttributeError(
            """ There is no attribute named rpacket_tracker in the simulation object passed. Try enabling the
                rpacket tracking in the configuration. To enable rpacket tracking see: https://tardis-sn.github.io/tardis/io/output/rpacket_tracking.html#How-to-Setup-the-Tracking-for-the-RPackets?"""
        )

    def generate_plot(self, theme="light"):
        """
        Creates an animated plotly plot showing the Montecarlo packets' trajectories.

        Parameters
        ----------
        theme : str, optional
            theme for the plot, by default "light"

        Returns
        -------
        plotly.graph_objs._figure.Figure
            plot containing the packets, photosphere and the shells.
        """
        self.fig = go.Figure()

        # getting velocity of different shells
        v_shells = self.sim.simulation_state.velocity.to_value(u.km / u.s)

        # getting coordinates and interactions of all packets
        (
            rpacket_x,
            rpacket_y,
            rpacket_interactions,
        ) = self.get_coordinates_multiple_packets(
            self.sim.transport.transport_state.rpacket_tracker_df.loc[
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
                marker=dict(color="#3366FF"),
            )
        )
        self.fig.add_trace(
            go.Scatter(
                **legend_props,
                name="Line Interaction",
                marker=dict(color="#FF3300"),
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
        r_track,
        mu_track,
        time,
        last_interaction_type,
        theta_initial=0,
    ):
        """
        Generates the coordinates of a single packet for its entire trajectory using the `nu` and `r` attributes.

        Parameters
        ----------
        r_track : pandas.core.series.Series
            radius of packet at different steps
        mu_track : pandas.core.series.Series
            mu or the cosine of radial angle of the packet at different steps
        time : astropy.units.quantity.Quantity
            time since the occurence of explosion
        last_interaction_type : pandas.core.series.Series
            last interaction type of the packet at different steps
        theta_initial : float, optional
            intial launch angle of packet from x axis as the packet starts from the photosphere, by default 0

        Returns
        -------
        list
            x coordinates of the packet at different steps
        list
            y coordinates of the packet at different steps
        list
            types of interactions occuring at different points
        """
        theta, rpacket_interactions = [], []

        # getting thetas at different steps of the packet movement

        for step_no in range(len(r_track)):
            # for the first step the packet is at photosphere, so theta will be equal to the intial angle we are launching the packet from
            if step_no == 0:
                theta.append(theta_initial)
            # for further steps we calculate thetas with the formula derived in the documentation
            else:
                curr_r = r_track[step_no]
                prev_r = r_track[step_no - 1]
                prev_mu = mu_track[step_no - 1]

                acos_mu = math.acos(prev_mu)
                sin_term = prev_r * math.sin(acos_mu) / curr_r
                new_theta = theta[-1] + acos_mu

                if curr_r < prev_r:
                    new_theta = new_theta - math.pi + math.asin(sin_term)
                else:
                    new_theta += math.asin(-1 * sin_term)
                theta.append(new_theta)

        # converting the thetas into x and y coordinates using radius as radius*cos(theta) and radius*sin(theta) respectively
        rpacket_x = (np.array(r_track)) * np.cos(np.array(theta)) * 1e-5 / time
        rpacket_y = (np.array(r_track)) * np.sin(np.array(theta)) * 1e-5 / time

        # adding interactions at different steps
        # using the change of slope of the trajectory line at different steps, we determine if an interactions happened or not.

        for step_no in range(len(r_track)):
            # when packet is at its starting and ending point in its trajectory, we consider it as no interaction
            if step_no == 0 or step_no == len(r_track) - 1:
                rpacket_interactions.append(0)
            else:
                # current slope is the slope of line from previous position of the packet to the current position
                rpacket_interactions.append(last_interaction_type[step_no])

        return rpacket_x, rpacket_y, rpacket_interactions

    def get_coordinates_multiple_packets(self, r_packet_tracker):
        """
        Generates an array of array containing x and y coordinates of multiple packets

        Parameters
        ----------
        r_packet_tracker : pandas.core.frame.DataFrame
            contains the rpacket_tracker_df dataframe with data of only specified no_of_packets

        Returns
        -------
        numpy.ndarray
            array of array containing x coordinates, y coordinates and the interactions for multiple packets
        """
        # for plotting packets at equal intervals throught the circle, we choose thetas distributed uniformly
        thetas = np.linspace(0, 2 * math.pi, self.no_of_packets + 1)
        all_rpackets_x_coords = []
        all_rpackets_y_coords = []
        all_rpackets_interactions_coords = []

        # getting coordinates and interaction arrays for all packets
        for packet_no in range(self.no_of_packets):
            (
                rpacket_x,
                rpacket_y,
                rpacket_interactions,
            ) = self.get_coordinates_with_theta_init(
                r_packet_tracker.loc[packet_no]["r"],
                r_packet_tracker.loc[packet_no]["mu"],
                self.sim.simulation_state.time_explosion.value,
                r_packet_tracker.loc[packet_no]["interaction_type"],
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

    def get_equal_array_size(self, rpacket_x, rpacket_y, interactions):
        """
        creates the coordinate arrays of different packets of same size. This is done for generating frames in animation.

        Parameters
        ----------
        rpacket_x : numpy.ndarray
            x coordinates of packets
        rpacket_y : numpy.ndarray
            y coordinates of packets
        interactions : numpy.ndarray
            interaction types of packets

        Returns
        -------
        numpy.ndarray
            normalized x coordinate array
        numpy.ndarray
            normalized y coordinate array
        numpy.ndarray
            normalized interaction types array
        int
            size of the biggest array among different packets

        """
        rpacket_step_no_array_max_size = max(list(map(len, rpacket_x)))

        for packet_no in range(len(rpacket_x)):
            padding_size = rpacket_step_no_array_max_size - len(
                rpacket_x[packet_no]
            )
            if padding_size > 0:
                x_padding = np.ones(padding_size) * rpacket_x[packet_no][-1]
                y_padding = np.ones(padding_size) * rpacket_y[packet_no][-1]
                interactions_padding = (
                    np.ones(padding_size) * interactions[packet_no][-1]
                )

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

    def get_frames(self, frame, rpacket_x, rpacket_y, interactions, theme):
        """
        Creates individual frames containing the go.Scatter objects for the animation.

        Parameters
        ----------
        frame : int
            current frame number
        rpacket_x : numpy.ndarray
            x coordinates array
        rpacket_y : numpy.ndarray
            y coordinates array
        interactions : numpy.ndarray
            interactions array
        theme : str
            theme for the plot

        Returns
        -------
        list
            list of go.Scatter objects for a particular frame number.
        """
        return [
            self.create_packet_scatter(
                packet_no, rpacket_x, rpacket_y, interactions, theme, frame=frame
            )
            for packet_no in range(len(rpacket_x))
        ]

    def create_packet_scatter(
        self, packet_no, rpacket_x, rpacket_y, interactions, theme, frame=None
    ):
        """
        Creates a Scatter object for a given packet with optional frame slicing.
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
            self.interaction_from_num[int(interaction)]["text"]
            for interaction in interaction_vals
        ]
        opacity_values = [
            self.interaction_from_num[int(interaction)]["opacity"]
            for interaction in interaction_vals
        ]
        color_values = [
            self.interaction_from_num[int(interaction)]["color"]
            for interaction in interaction_vals
        ]

        return go.Scatter(
            x=x_vals,
            y=y_vals,
            mode="markers+lines",
            name=f"Packet {str(packet_no + 1)}",
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

    def get_slider_steps(self, rpacket_max_array_size):
        """
        Generates different steps in the timeline slider for different frames in the animated plot.

        Parameters
        ----------
        rpacket_max_array_size : int
            maximum size of coordinate array among all the packets.

        Returns
        -------
        list
            list of dictionaries of different steps for different frames.
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
