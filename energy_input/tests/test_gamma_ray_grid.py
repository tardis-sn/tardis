import pytest
import numpy.testing as npt
import numpy as np

from tardis.energy_input.gamma_ray_grid import (
    calculate_distance_radial,
    distance_trace,
    move_packet,
)


@pytest.mark.xfail(reason="To be implemented")
def test_calculate_distance_radial():
    """Test the radial distance calculation"""
    assert False


@pytest.mark.xfail(reason="To be implemented")
def test_distance_trace():
    """Test the distance trace"""
    assert False


def test_move_packet(basic_gamma_ray):
    """

    Parameters
    ----------
    basic_gamma_ray : GammaRay object
    """
    packet = basic_gamma_ray
    distance = 1.0e15

    new = packet.location + distance * packet.direction

    actual = move_packet(packet, distance)

    npt.assert_almost_equal(actual.location, new)
