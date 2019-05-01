from numba import int64, float64
from numba import jitclass, njit


rpacket_spec = [
    ('r', float64),
    ('mu', float64),
    ('tau_event', float64)
    ('nu_line', float64)
    ('d_line', float64)  # Distance to electron event. 
    ('d_electron', float64) #/**< Distance to line event. */
    ('d_boundary', float64) # distance to boundary 
    ('shell_id', int64)
]

@jitclass(rpacket_rspec)
class RPacket(object):
    def __init__(self, r, mu, tau_event):
        self.r = r
        self.mu = mu
        self.nu = nu
        self.tau_event = tau_event
        self.nu_line = nu_line

    @staticmethod
    def get_tau_event():
        return -np.log(np.random.random())