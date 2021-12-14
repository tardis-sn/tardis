from enum import IntEnum
import numpy as np


class GXPhotonStatus(IntEnum):
    BETA_DECAY = -1
    COMPTON_SCATTER = 0
    PHOTOABSORPTION = 1
    PAIR_CREATION = 2
    IN_PROCESS = 3
    END = 4


class GXPhoton(object):
    """
    Gamma ray or X ray object with location, direction, energy, time and optical depth

    Attributes
    ----------
    location : SphericalVector object
             GXPhoton position vector
    direction : SphericalVector object
             GXPhoton direction vector (unitary)
    energy : float64
             GXPhoton energy
    status : InteractionType
             GXPhoton status
    shell : int64
             GXPhoton shell location index
    """

    def __init__(self, location, direction, energy, status, shell, activity):
        self.location = location
        self.direction = direction
        self.energy = energy
        self.status = status
        self.shell = shell
        self.time_created = 0
        self.time_current = 0
        self.tau = -np.log(np.random.random())
        self.activity = activity
