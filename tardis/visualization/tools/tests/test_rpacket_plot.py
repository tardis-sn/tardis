"""Tests for RPacketPlotter Plots"""

import math
import astropy.units as u
import pytest
import numpy as np
import numpy.testing as npt
from tardis.visualization import RPacketPlotter


@pytest.mark.rpacket_tracking
class TestRPacketPlotter:
    """Test the RPacketPlotter class."""

    def test_get_coordinates_with_theta_init(self, simulation_rpacket_tracking):
        """
        Test for the get_coordinates_with_theta_init method.

        Parameters
        ----------
        simulation_rpacket_tracking : tardis.simulation.base.Simulation
            Simulation object.
        """
        rpacket_plotter = RPacketPlotter.from_simulation(
            simulation_rpacket_tracking
        )
        single_packet_df = simulation_rpacket_tracking.transport.transport_state.rpacket_tracker_df.loc[
            0
        ]
        (
            rpacket_x,
            rpacket_y,
            rpacket_interactions,
        ) = rpacket_plotter.get_coordinates_with_theta_init(
            single_packet_df["r"],
            single_packet_df["mu"],
            simulation_rpacket_tracking.simulation_state.time_explosion.value,
            single_packet_df["interaction_type"],
            0,
        )

        # checking the length of coordinate array
        assert len(rpacket_x) == single_packet_df.shape[0]
        assert len(rpacket_y) == single_packet_df.shape[0]
        assert len(rpacket_interactions) == single_packet_df.shape[0]

        # compairing the radius obtained from x and y coordinates using the formula r^2 = x^2 + y^2,  with the radius present in the rpacket_tracker
        radius_array = np.sqrt(rpacket_x**2 + rpacket_y**2)
        expected_radius_array = (
            (single_packet_df["r"].to_numpy())
            * 1e-5
            / simulation_rpacket_tracking.simulation_state.time_explosion.value
        )
        npt.assert_allclose(radius_array, expected_radius_array)

    @pytest.mark.parametrize("no_of_packets", [2, 5, 10])
    def test_get_coordinates_multiple_packets(
        self, simulation_rpacket_tracking, no_of_packets
    ):
        """
        Test for the get_coordinates_multiple_packets method.

        Parameters
        ----------
        simulation_rpacket_tracking : tardis.simulation.base.Simulation
            Simulation object.
        no_of_packets : int
            Number of Packets.
        """
        rpacket_plotter = RPacketPlotter.from_simulation(
            simulation_rpacket_tracking, no_of_packets=no_of_packets
        )
        multiple_packet_df = simulation_rpacket_tracking.transport.transport_state.rpacket_tracker_df.loc[
            0 : (no_of_packets - 1)
        ]
        (
            rpackets_x,
            rpackets_y,
            rpackets_interactions,
        ) = rpacket_plotter.get_coordinates_multiple_packets(multiple_packet_df)
        # checking the length of coordiunate array
        assert len(rpackets_x) == no_of_packets
        assert len(rpackets_y) == no_of_packets
        assert len(rpackets_interactions) == no_of_packets

        thetas = np.linspace(0, 2 * math.pi, no_of_packets + 1)

        # checking coordinates of every packet
        for rpacket in range(no_of_packets):
            single_packet_df = simulation_rpacket_tracking.transport.transport_state.rpacket_tracker_df.loc[
                rpacket
            ]
            (
                expected_rpacket_x,
                expected_rpacket_y,
                expected_rpacket_interactions,
            ) = rpacket_plotter.get_coordinates_with_theta_init(
                single_packet_df["r"],
                single_packet_df["mu"],
                simulation_rpacket_tracking.simulation_state.time_explosion.value,
                single_packet_df["interaction_type"],
                thetas[rpacket],
            )
            npt.assert_array_equal(rpackets_x[rpacket], expected_rpacket_x)
            npt.assert_array_equal(rpackets_y[rpacket], expected_rpacket_y)
            npt.assert_array_equal(
                rpackets_interactions[rpacket], expected_rpacket_interactions
            )

    @pytest.mark.parametrize("no_of_packets", [2, 5, 10])
    def test_get_equal_array_size(
        self, simulation_rpacket_tracking, no_of_packets
    ):
        """
        Test for the get_equal_array_size method.

        Parameters
        ----------
        simulation_rpacket_tracking: tardis.simulation.base.Simulation
            Simulation object.
        no_of_packets : int
            Number of Packets
        """
        rpacket_plotter = RPacketPlotter.from_simulation(
            simulation_rpacket_tracking, no_of_packets=no_of_packets
        )
        multiple_packet_df = simulation_rpacket_tracking.transport.transport_state.rpacket_tracker_df.loc[
            0 : (no_of_packets - 1)
        ]
        (
            multiple_packet_x,
            multiple_packet_y,
            multiple_packet_interaction,
        ) = rpacket_plotter.get_coordinates_multiple_packets(multiple_packet_df)
        (
            rpackets_x,
            rpackets_y,
            rpackets_interactions,
            max_array_size,
        ) = rpacket_plotter.get_equal_array_size(
            multiple_packet_x, multiple_packet_y, multiple_packet_interaction
        )

        expected_max_array_size = max(list(map(len, multiple_packet_x)))

        # checking the calculated max_array_size value
        assert expected_max_array_size == max_array_size

        # checking the final array for all packets
        for rpacket in range(no_of_packets):
            expected_rpacket_x = np.append(
                multiple_packet_x[rpacket],
                multiple_packet_x[rpacket][-1]
                * np.ones(
                    [expected_max_array_size - len(multiple_packet_x[rpacket])]
                ),
            )
            expected_rpacket_y = np.append(
                multiple_packet_y[rpacket],
                multiple_packet_y[rpacket][-1]
                * np.ones(
                    [expected_max_array_size - len(multiple_packet_y[rpacket])]
                ),
            )
            expected_rpacket_interactions = np.append(
                multiple_packet_interaction[rpacket],
                multiple_packet_interaction[rpacket][-1]
                * np.ones(
                    [
                        expected_max_array_size
                        - len(multiple_packet_interaction[rpacket])
                    ]
                ),
            )

            npt.assert_array_equal(expected_rpacket_x, rpackets_x[rpacket])
            npt.assert_array_equal(expected_rpacket_y, rpackets_y[rpacket])
            npt.assert_array_equal(
                expected_rpacket_interactions, rpackets_interactions[rpacket]
            )

    @pytest.mark.parametrize("no_of_packets", [2, 5, 10])
    @pytest.mark.parametrize("theme", ["light", "dark"])
    def test_get_frames(
        self, simulation_rpacket_tracking, no_of_packets, theme
    ):
        """
        Test for the get_frames method.

        Parameters
        ----------
        simulation_rpacket_tracking : tardis.simulation.base.Simulation
            Simulation object.
        no_of_packets : int
            Number of Packets
        theme : str
            Theme of plot.
        """
        rpacket_plotter = RPacketPlotter.from_simulation(
            simulation_rpacket_tracking, no_of_packets=no_of_packets
        )
        multiple_packet_df = simulation_rpacket_tracking.transport.transport_state.rpacket_tracker_df.loc[
            0 : (no_of_packets - 1)
        ]
        (
            multiple_packet_x,
            multiple_packet_y,
            multiple_packet_interaction,
        ) = rpacket_plotter.get_coordinates_multiple_packets(multiple_packet_df)
        (
            rpackets_x,
            rpackets_y,
            rpackets_interactions,
            rpacket_array_max_size,
        ) = rpacket_plotter.get_equal_array_size(
            multiple_packet_x, multiple_packet_y, multiple_packet_interaction
        )

        # checking data inside all frames
        for frame in range(rpacket_array_max_size + 1):
            frames = rpacket_plotter.get_frames(
                frame, rpackets_x, rpackets_y, rpackets_interactions, theme
            )
            assert len(frames) == no_of_packets
            for index, packet_frame in enumerate(frames):
                expected_x = rpackets_x[index].tolist()[0:frame]
                expected_y = rpackets_y[index].tolist()[0:frame]
                npt.assert_allclose(expected_x, packet_frame.x)
                npt.assert_allclose(expected_y, packet_frame.y)

    @pytest.mark.parametrize("max_step_size", [10, 30, 50])
    def test_get_slider_steps(self, simulation_rpacket_tracking, max_step_size):
        rpacket_plotter = RPacketPlotter.from_simulation(
            simulation_rpacket_tracking
        )
        slider_steps = rpacket_plotter.get_slider_steps(max_step_size)
        for index, step in enumerate(slider_steps):
            assert step["args"][0][0] == index
            assert step["label"] == index

    @pytest.mark.parametrize("no_of_packets", [2, 5, 10])
    @pytest.mark.parametrize("theme", ["light", "dark"])
    def test_generate_plot(
        self, simulation_rpacket_tracking, no_of_packets, theme
    ):
        """
        Test for the generate_plot method.

        Parameters
        ----------
        simulation_rpacket_tracking : tardis.simulation.base.Simulation
            Simulation object.
        no_of_packets : int
            Number of Packets
        theme : str
            Theme of plot.
        """
        rpacket_plotter = RPacketPlotter.from_simulation(
            simulation_rpacket_tracking, no_of_packets=no_of_packets
        )
        multiple_packet_df = simulation_rpacket_tracking.transport.transport_state.rpacket_tracker_df.loc[
            0 : (no_of_packets - 1)
        ]
        (
            multiple_packet_x,
            multiple_packet_y,
            multiple_packet_interaction,
        ) = rpacket_plotter.get_coordinates_multiple_packets(multiple_packet_df)
        (
            rpackets_x,
            rpackets_y,
            rpackets_interactions,
            rpacket_array_max_size,
        ) = rpacket_plotter.get_equal_array_size(
            multiple_packet_x, multiple_packet_y, multiple_packet_interaction
        )

        fig = rpacket_plotter.generate_plot(theme=theme)

        shell_radii = (
            simulation_rpacket_tracking.simulation_state.velocity.to_value(
                u.km / u.s
            )
        )

        # testing the shells and the photosphere
        for index, shell in enumerate(fig.layout.shapes):
            assert shell.type == "circle"
            assert shell.x0 == -1 * shell_radii[index]
            assert shell.x1 == shell_radii[index]
            assert shell.y0 == -1 * shell_radii[index]
            assert shell.y1 == shell_radii[index]

        # testing the non-animated static plot
        for packet_no in range(no_of_packets):
            packet = fig.data[packet_no]
            assert packet.name == "Packet " + str(packet_no + 1)
            npt.assert_allclose(packet.x, rpackets_x[packet_no])
            npt.assert_allclose(packet.y, rpackets_y[packet_no])
            assert list(packet.marker.color) == [
                rpacket_plotter.interaction_from_num[
                    int(rpackets_interactions[packet_no][step_no])
                ]["color"]
                for step_no in range(len(rpackets_x[packet_no]))
            ]

        # testing the frames generated in the animation
        for frame in fig.frames:
            for index, frame in enumerate(fig.frames):
                expected_frame = rpacket_plotter.get_frames(
                    index, rpackets_x, rpackets_y, rpackets_interactions, theme
                )
                for packet in range(no_of_packets):
                    assert frame.data[packet] == expected_frame[packet]
