import pytest
import numpy.testing as npt
import numpy as np


from tardis.transport.r_packet_transport_3d import (
    move_packet_3d,
)


def test_move_packet_3d(basic_gamma_ray):
    """

    Parameters
    ----------
    basic_gamma_ray : GammaRay object
    """
    packet = basic_gamma_ray
    distance = 1.0e15

    new = packet.location + distance * packet.direction

    actual = move_packet_3d(packet, distance)

    npt.assert_almost_equal(actual.location, new)
