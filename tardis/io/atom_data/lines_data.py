import numpy as np
from astropy.units import Quantity


class LinesData:
    def __init__(self, lines):
        # Convert wavelengths to CGS
        lines["wavelength_cm"] = Quantity(
            lines["wavelength"], "angstrom"
        ).cgs.value
        self.data = lines

    def active_data(self):
        pass


class LineList:
    def __init__(self, linelist):
        self.data = linelist
