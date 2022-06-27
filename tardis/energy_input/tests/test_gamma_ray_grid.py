import pytest
import numpy.testing as npt
import numpy as np

from tardis.transport.geometry.calculate_distances_pseudo_3D import (
    calculate_distance_boundary_pseudo_3D,
    distance_trace_gamma,
)

from tardis.transport.r_packet_transport_pseudo_3D import (
    move_packet_pseudo_3D
)


@pytest.mark.xfail(reason="To be implemented")
def test_calculate_distance_boundary_pseudo_3D():
    """Test the radial distance calculation"""
    assert False


@pytest.mark.xfail(reason="To be implemented")
def test_distance_trace_gamma():
    """Test the distance trace"""
    assert False


def test_move_packet_pseudo_3D(basic_gamma_ray):
    """

    Parameters
    ----------
    basic_gamma_ray : GammaRay object
    """
    packet = basic_gamma_ray
    distance = 1.0e15

    new = packet.location + distance * packet.direction

    actual = move_packet_pseudo_3D(packet, distance)

    npt.assert_almost_equal(actual.location, new)
