import numpy.testing as npt

from tardis.energy_input.transport.gamma_ray_grid import (
    move_packet,
)


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
