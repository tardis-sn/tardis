import astropy.units as u
import plotly.express as px
import plotly.graph_objects as go
import math
import pandas as pd
import numpy as np
import random


class RPacketPlotter:
    def __init__(self, sim, no_of_packets):
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
        self.fig = go.Figure()
        v_shells = self.sim.model.velocity.to_value(u.km / u.s)
        xs, ys, ints = self.get_coordinates_multiple_packets(
            self.sim.runner.rpacket_tracker_df.loc[0 : (self.no_of_packets)],
        )
        xs, ys, ints, max_size = self.get_equal_array_size(xs, ys, ints)
        # Set axes properties
        self.fig.update_xaxes(
            range=[-1.1 * v_shells[-1], 1.1 * v_shells[-1]],
            title="VELOCITY (KM/S)",
            exponentformat="none",
            linecolor="#555",
            gridcolor="#fafafa",
            zerolinecolor="#fafafa",
        )
        self.fig.update_yaxes(
            range=[-1.1 * v_shells[-1], 1.1 * v_shells[-1]],
            title="VELOCITY (KM/S)",
            exponentformat="none",
            linecolor="#555",
            gridcolor="#fafafa",
            zerolinecolor="#fafafa",
        )

        # adding the shells and photosphere
        shell_shapes = {}
        for i in range(len(self.sim.model.radius.value)):
            if i == 0:
                # photosphere
                self.fig.add_shape(
                    type="circle",
                    xref="x",
                    yref="y",
                    x0=-1 * v_shells[i],
                    y0=-1 * v_shells[i],
                    x1=v_shells[i],
                    y1=v_shells[i],
                    line_color="black",
                    fillcolor="darkgrey",
                    opacity=1,
                )
            elif i == (len(self.sim.model.radius.value) - 1):
                # outermost shell
                self.fig.add_shape(
                    type="circle",
                    xref="x",
                    yref="y",
                    x0=-1 * v_shells[i],
                    y0=-1 * v_shells[i],
                    x1=v_shells[i],
                    y1=v_shells[i],
                    line_color="black",
                    opacity=1,
                )
            else:
                # remaining shells
                self.fig.add_shape(
                    type="circle",
                    xref="x",
                    yref="y",
                    x0=-1 * v_shells[i],
                    y0=-1 * v_shells[i],
                    x1=v_shells[i],
                    y1=v_shells[i],
                    line_color="black",
                    opacity=0.1,
                )

        # Adding packet trajectory

        for i in range(len(xs)):
            self.fig.add_trace(
                go.Scatter(
                    x=xs[i],
                    y=ys[i],
                    mode="markers+lines",
                    name="Packet " + str(i + 1),
                    showlegend=False,
                    hovertemplate="<b>X</b>: %{x}"
                    + "<br><b>Y</b>: %{y}<br>"
                    + "<b>Last Interaction: %{text}</b>",
                    text=[
                        self.interaction_from_num.get(ints[i][j])
                        for j in range(len(xs[i]))
                    ],
                    line=dict(color="darkslategrey"),
                    marker=dict(
                        opacity=[
                            self.interaction_opacity_from_num.get(ints[i][j])
                            for j in range(len(xs[i]))
                        ],
                        color=[
                            self.interaction_color_from_num.get(ints[i][j])
                            for j in range(len(xs[i]))
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
            height=900,
            title="Packet Trajectories",
            title_font_color="#444",
            updatemenus=[
                dict(
                    type="buttons",
                    pad=dict(t=750),
                    buttons=[dict(label="Play", method="animate", args=[None])],
                )
            ],
        )

        self.fig.frames = [
            go.Frame(data=self.get_frames(frame, xs, ys, ints))
            for frame in range(max_size + 1)
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
        xs, ys, theta, ints = [], [], [], []

        for i in range(len(r_track)):
            if i == 0:
                theta.append(theta_initial)
            else:
                if r_track[i] < r_track[i - 1]:
                    theta.append(
                        theta[-1]
                        - math.pi
                        + math.asin(
                            r_track[i - 1]
                            * math.sin(math.acos(mu_track[i - 1]))
                            / r_track[i]
                        )
                        + math.acos(mu_track[i - 1])
                    )
                else:
                    theta.append(
                        theta[-1]
                        + math.asin(
                            -1
                            * r_track[i - 1]
                            * math.sin(math.acos(mu_track[i - 1]))
                            / r_track[i]
                        )
                        + math.acos(mu_track[i - 1])
                    )

        xs = (np.array(r_track)) * np.cos(np.array(theta)) * 1e-5 / time
        ys = (np.array(r_track)) * np.sin(np.array(theta)) * 1e-5 / time

        for i in range(len(r_track)):
            if i == 0 or i == len(r_track) - 1:
                ints.append(0)
            else:
                s0 = (ys[i] - ys[i - 1]) / (xs[i] - xs[i - 1])
                s1 = (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i])
                if math.isclose(s0, s1, rel_tol=1e-11):
                    ints.append(0)
                else:
                    ints.append(last_interaction_type[i])

        return xs, ys, ints

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
        x = []
        y = []
        inters = []
        for i in range(self.no_of_packets):
            xs, ys, ints = self.get_coordinates_with_theta_init(
                r_packet_tracker.loc[i]["r"],
                r_packet_tracker.loc[i]["mu"],
                self.sim.model.time_explosion.value,
                r_packet_tracker.loc[i]["interaction_type"],
                thetas[i],
            )
            x.append(xs)
            y.append(ys)
            inters.append(ints)
        return (
            np.array(x, dtype="object"),
            np.array(y, dtype="object"),
            np.array(inters, dtype="object"),
        )

    def get_equal_array_size(self, xs, ys, interactions):
        """
        creates the coordinate arrays of different packets of same size. This is done for generating frames in animation.

        Parameters
        ----------
        xs : numpy.ndarray
            x coordinates of packets
        ys : numpy.ndarray
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
        max_size = max(list(map(len, xs)))
        for i in range(len(xs)):
            xs[i] = np.append(
                xs[i], xs[i][-1] * np.ones([max_size - len(xs[i])])
            )
            ys[i] = np.append(
                ys[i], ys[i][-1] * np.ones([max_size - len(ys[i])])
            )
            interactions[i] = np.append(
                interactions[i],
                interactions[i][-1]
                * np.ones([max_size - len(interactions[i])]),
            )
        return xs, ys, interactions, max_size

    # creating frames for animation
    def get_frames(self, frame, xs, ys, interactions):
        """
        Creates individual frames containing the go.Scatter objects for the animation.

        Parameters
        ----------
        frame : int
            current frame number
        xs : numpy.ndarray
            x coordinates array
        ys : numpy.ndarray
            y coordinates array
        interactions : numpy.ndarray
            interactions array

        Returns
        -------
        list
            list of go.Scatter objects for a particular frame number.
        """
        frames = []
        for i in range(len(xs)):
            frames.append(
                go.Scatter(
                    x=xs[i].tolist()[0:frame],
                    y=ys[i].tolist()[0:frame],
                    mode="markers+lines",
                    name="Packet " + str(i + 1),
                    showlegend=False,
                    hovertemplate="<b>X</b>: %{x}"
                    + "<br><b>Y</b>: %{y}<br>"
                    + "<b>Last Interaction: %{text}</b>",
                    text=[
                        self.interaction_from_num.get(interactions[i][j])
                        for j in range(len(xs[i]))
                    ],
                    line=dict(color="darkslategrey"),
                    marker=dict(
                        opacity=[
                            self.interaction_opacity_from_num.get(
                                interactions[i][j]
                            )
                            for j in range(len(xs[i]))
                        ],
                        color=[
                            self.interaction_color_from_num.get(
                                interactions[i][j]
                            )
                            for j in range(len(xs[i]))
                        ],
                    ),
                )
            )
        return frames
