from numba import float64, jitclass
from tardis import constants as const


C_SPEED_OF_LIGHT = const.c.to('cm/s').value

numba_model_spec = [
    ('r_inner', float64[:]),
    ('r_outer', float64[:]),
    ('time_explosion', float64),
    ('electron_density', float64[:]),
    ('ct', float64),
]


@jitclass(numba_model_spec)
class NumbaModel(object):

    def __init__(self, r_inner, r_outer, time_explosion, electron_density):
        """
        Model for the Numba mode

        Parameters
        ----------
        r_inner: numpy.ndarray
        r_outer: numpy.ndarray
        time_explosion: float
        electron_density: numpy.ndarray
        """
        self.r_inner = r_inner
        self.r_outer = r_outer
        self.time_explosion = time_explosion
        self.electron_density = electron_density


packet_collection_spec = [
    ('packets_nu', float64[:]),
    ('packets_mu', float64[:]),
    ('packets_energy', float64[:]),
    ('packets_output_nu', float64[:]),
    ('packets_output_energy', float64[:]),
]

@jitclass(packet_collection_spec)
class PacketCollection(object):
    def __init__(self, packets_input_nu, packets_input_mu, packets_energy,
                 packets_output_nu, packets_output_energy):
        self.packets_input_nu = packets_input_nu
        self.packets_input_mu = packets_input_mu
        self.packets_energy = packets_energy
        self.packets_output_nu = packets_output_nu
        self.packets_output_energy = packets_output_energy




class Estimators(object):
    pass

