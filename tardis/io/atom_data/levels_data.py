import numpy as np
from astropy.units import Quantity


class LevelsData:
    def __init__(self, levels):
        # Convert energy to CGS
        levels.loc[:, "energy"] = Quantity(
            levels["energy"].values, "eV"
        ).cgs.value
        self.data = levels
        # Cast to float so that Numba can use the values in numpy functions
        self.data.energy = self.data.energy.astype(np.float64)

    def active_data(self):
        pass
