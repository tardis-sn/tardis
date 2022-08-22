import astropy.units as u
import plotly.express as px
import plotly.graph_objects as go
import math
import pandas as pd
import numpy as np
import random


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
        self.interaction_from_num = {
            0: "No Interaction",
            1: "EScattering",
            2: "Line",
        }
        self.interaction_color_from_num = {
            0: "darkslategrey",
            1: "#3366FF",
            2: "#FF3300",
        }
        self.interaction_opacity_from_num = {0: 0, 1: 1, 2: 1}

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
        SDECPlotter
        """
        if hasattr(sim.runner, "rpacket_tracker_df"):
            if len(sim.runner.rpacket_tracker) >= no_of_packets:
                return cls(sim, no_of_packets)
            else:
                raise ValueError(
                    """
                    no_of_packets specified are more than the actual no of packets in the model.
                    """
                )
        else:
            raise NameError(
                """
                There is no attribute named rpacket_tracker. Try enabling the 
                rpacket tracking in the configuration.
                """
            )

    def generate_plot(self):
        """
        Creates an animated plotly plot showing the Montecarlo packets' trajectories.

        Returns
        -------
        plotly.graph_objs._figure.Figure
            plot containing the packets, photosphere and the shells.
        """

        self.fig = go.Figure()
        v_shells = self.sim.model.velocity.to_value(u.km / u.s)
        (
            rpacket_x,
            rpacket_y,
            rpacket_interactions,
        ) = self.get_coordinates_multiple_packets(
            self.sim.runner.rpacket_tracker_df.loc[0 : (self.no_of_packets)],
        )
        (
            rpacket_x,
            rpacket_y,
            rpacket_interactions,
            rpacket_array_max_size,
        ) = self.get_equal_array_size(
            rpacket_x, rpacket_y, rpacket_interactions
        )
        # Set axes properties
        self.fig.update_xaxes(
            scaleanchor="y",
            scaleratio=1,
            range=[-1.1 * v_shells[-1], 1.1 * v_shells[-1]],
            title="Velocity (km/s)",
            exponentformat="none",
            linecolor="#555",
            gridcolor="#fafafa",
            zerolinecolor="#fafafa",
        )
        self.fig.update_yaxes(
            range=[-1.1 * v_shells[-1], 1.1 * v_shells[-1]],
            title="Velocity (km/s)",
            exponentformat="none",
            linecolor="#555",
            gridcolor="#fafafa",
            zerolinecolor="#fafafa",
        )

        # adding the shells and photosphere
        for shell_no in range(len(self.sim.model.radius.value)):
            if shell_no == 0:
                # photosphere
                self.fig.add_shape(
                    type="circle",
                    xref="x",
                    yref="y",
                    x0=-1 * v_shells[shell_no],
                    y0=-1 * v_shells[shell_no],
                    x1=v_shells[shell_no],
                    y1=v_shells[shell_no],
                    line_color="black",
                    fillcolor="darkgrey",
                    opacity=1,
                )
            elif shell_no == (len(self.sim.model.radius.value) - 1):
                # outermost shell
                self.fig.add_shape(
                    type="circle",
                    xref="x",
                    yref="y",
                    x0=-1 * v_shells[shell_no],
                    y0=-1 * v_shells[shell_no],
                    x1=v_shells[shell_no],
                    y1=v_shells[shell_no],
                    line_color="black",
                    opacity=1,
                )
            else:
                # remaining shells
                self.fig.add_shape(
                    type="circle",
                    xref="x",
                    yref="y",
                    x0=-1 * v_shells[shell_no],
                    y0=-1 * v_shells[shell_no],
                    x1=v_shells[shell_no],
                    y1=v_shells[shell_no],
                    line_color="black",
                    opacity=0.1,
                )

        # Adding packet trajectory

        for packet_no in range(len(rpacket_x)):
            self.fig.add_trace(
                go.Scatter(
                    x=rpacket_x[packet_no],
                    y=rpacket_y[packet_no],
                    mode="markers+lines",
                    name="Packet " + str(packet_no + 1),
                    showlegend=False,
                    hovertemplate="<b>X</b>: %{x}"
                    + "<br><b>Y</b>: %{y}<br>"
                    + "<b>Last Interaction: %{text}</b>",
                    text=[
                        self.interaction_from_num.get(
                            rpacket_interactions[packet_no][step_no]
                        )
                        for step_no in range(len(rpacket_x[packet_no]))
                    ],
                    line=dict(color="darkslategrey"),
                    marker=dict(
                        opacity=[
                            self.interaction_opacity_from_num.get(
                                rpacket_interactions[packet_no][step_no]
                            )
                            for step_no in range(len(rpacket_x[packet_no]))
                        ],
                        color=[
                            self.interaction_color_from_num.get(
                                rpacket_interactions[packet_no][step_no]
                            )
                            for step_no in range(len(rpacket_x[packet_no]))
                        ],
                    ),
                )
            )

        # adding legends
        self.fig.add_trace(
            go.Scatter(
                x=[9999999],
                y=[0],
                legendgroup="a",
                opacity=1,
                legendgrouptitle=dict(
                    font=dict(color="#444"), text="Interaction Type:"
                ),
                mode="lines+markers",
                name="Escattering",
                hoverlabel=dict(font=dict(color="#222")),
                marker=dict(color="#3366FF"),
            )
        )
        self.fig.add_trace(
            go.Scatter(
                x=[9999999],
                y=[0],
                legendgroup="a",
                opacity=1,
                mode="lines+markers",
                name="Line",
                marker=dict(color="#FF3300"),
            )
        )

        # Set figure size
        self.fig.layout.plot_bgcolor = "#fafafa"
        self.fig.layout.paper_bgcolor = "#fafafa"

        self.fig.update_layout(
            width=900,
            height=700,
            title="Packet Trajectories",
            title_font_color="#444",
            updatemenus=[
                dict(
                    type="buttons",
                    y=-0.1,
                    buttons=[dict(label="Play", method="animate", args=[None])],
                )
            ],
        )

        self.fig.frames = [
            go.Frame(
                data=self.get_frames(
                    frame, rpacket_x, rpacket_y, rpacket_interactions
                )
            )
            for frame in range(rpacket_array_max_size + 1)
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
        rpacket_x, rpacket_y, theta, rpacket_interactions = [], [], [], []

        for step_no in range(len(r_track)):
            if step_no == 0:
                theta.append(theta_initial)
            else:
                if r_track[step_no] < r_track[step_no - 1]:
                    theta.append(
                        theta[-1]
                        - math.pi
                        + math.asin(
                            r_track[step_no - 1]
                            * math.sin(math.acos(mu_track[step_no - 1]))
                            / r_track[step_no]
                        )
                        + math.acos(mu_track[step_no - 1])
                    )
                else:
                    theta.append(
                        theta[-1]
                        + math.asin(
                            -1
                            * r_track[step_no - 1]
                            * math.sin(math.acos(mu_track[step_no - 1]))
                            / r_track[step_no]
                        )
                        + math.acos(mu_track[step_no - 1])
                    )

        rpacket_x = (np.array(r_track)) * np.cos(np.array(theta)) * 1e-5 / time
        rpacket_y = (np.array(r_track)) * np.sin(np.array(theta)) * 1e-5 / time

        for step_no in range(len(r_track)):
            if step_no == 0 or step_no == len(r_track) - 1:
                rpacket_interactions.append(0)
            else:
                current_slope = (
                    rpacket_y[step_no] - rpacket_y[step_no - 1]
                ) / (rpacket_x[step_no] - rpacket_x[step_no - 1])
                next_slope = (rpacket_y[step_no + 1] - rpacket_y[step_no]) / (
                    rpacket_x[step_no + 1] - rpacket_x[step_no]
                )
                if math.isclose(current_slope, next_slope, rel_tol=1e-11):
                    rpacket_interactions.append(0)
                else:
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

        thetas = np.linspace(0, 2 * math.pi, self.no_of_packets + 1)
        rpackets_x = []
        rpackets_y = []
        rpackets_interactions = []
        for packet_no in range(self.no_of_packets):
            (
                rpacket_x,
                rpacket_y,
                rpacket_interactions,
            ) = self.get_coordinates_with_theta_init(
                r_packet_tracker.loc[packet_no]["r"],
                r_packet_tracker.loc[packet_no]["mu"],
                self.sim.model.time_explosion.value,
                r_packet_tracker.loc[packet_no]["interaction_type"],
                thetas[packet_no],
            )
            rpackets_x.append(rpacket_x)
            rpackets_y.append(rpacket_y)
            rpackets_interactions.append(rpacket_interactions)
        return (
            np.array(rpackets_x, dtype="object"),
            np.array(rpackets_y, dtype="object"),
            np.array(rpackets_interactions, dtype="object"),
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
        rpacket_arraystep_nomax_size = max(list(map(len, rpacket_x)))
        for packet_no in range(len(rpacket_x)):
            rpacket_x[packet_no] = np.append(
                rpacket_x[packet_no],
                rpacket_x[packet_no][-1]
                * np.ones(
                    [rpacket_arraystep_nomax_size - len(rpacket_x[packet_no])]
                ),
            )
            rpacket_y[packet_no] = np.append(
                rpacket_y[packet_no],
                rpacket_y[packet_no][-1]
                * np.ones(
                    [rpacket_arraystep_nomax_size - len(rpacket_y[packet_no])]
                ),
            )
            interactions[packet_no] = np.append(
                interactions[packet_no],
                interactions[packet_no][-1]
                * np.ones(
                    [
                        rpacket_arraystep_nomax_size
                        - len(interactions[packet_no])
                    ]
                ),
            )
        return rpacket_x, rpacket_y, interactions, rpacket_arraystep_nomax_size

    # creating frames for animation
    def get_frames(self, frame, rpacket_x, rpacket_y, interactions):
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

        Returns
        -------
        list
            list of go.Scatter objects for a particular frame number.
        """
        frames = []
        for packet_no in range(len(rpacket_x)):
            frames.append(
                go.Scatter(
                    x=rpacket_x[packet_no].tolist()[0:frame],
                    y=rpacket_y[packet_no].tolist()[0:frame],
                    mode="markers+lines",
                    name="Packet " + str(packet_no + 1),
                    showlegend=False,
                    hovertemplate="<b>X</b>: %{x}"
                    + "<br><b>Y</b>: %{y}<br>"
                    + "<b>Last Interaction: %{text}</b>",
                    text=[
                        self.interaction_from_num.get(
                            interactions[packet_no][step_no]
                        )
                        for step_no in range(len(rpacket_x[packet_no]))
                    ],
                    line=dict(color="darkslategrey"),
                    marker=dict(
                        opacity=[
                            self.interaction_opacity_from_num.get(
                                interactions[packet_no][step_no]
                            )
                            for step_no in range(len(rpacket_x[packet_no]))
                        ],
                        color=[
                            self.interaction_color_from_num.get(
                                interactions[packet_no][step_no]
                            )
                            for step_no in range(len(rpacket_x[packet_no]))
                        ],
                    ),
                )
            )
        return frames
