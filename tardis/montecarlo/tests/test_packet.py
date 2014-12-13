import pytest

from tardis.montecarlo.tests import packet_test_wrapper

@pytest.fixture
def rpacket():
    rpacket_dict = {'close_line': 0,
     'current_shell_id': 1,
     'd_boundary': 1.0,
     'd_electron': 1.0,
     'd_line': 1.0,
     'energy': 1.0,
     'last_line': 0,
     'mu': 1.0,
     'next_line_id': 1,
     'next_shell_id': 1,
     'nu': 1.0,
     'nu_line': 1.0,
     'r': 1.0,
     'recently_crossed_boundary': 0,
     'tau_event': 1.0,
     'virtual_packet': 1,
     'virtual_packet_flag': 1}
    return packet_test_wrapper.RPacket(rpacket_dict)

def test_set_get_nu(rpacket):
    assert False