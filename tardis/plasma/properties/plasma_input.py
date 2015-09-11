import numpy as np
import pandas as pd
from tardis.plasma.properties.base import (Input, ArrayInput, DataFrameInput)

__all__ = ['TRadiative', 'DilutionFactor', 'AtomicData', 'Abundance', 'Density',
           'TimeExplosion', 'JBlues', 'LinkTRadTElectron']

class TRadiative(ArrayInput):
    """
    Outputs:
    t_rad : Numpy Array
    """
    outputs = ('t_rad',)
    latex_name = ('T_{\\textrm{rad}}',)

class DilutionFactor(ArrayInput):
    """
    Outputs:
    w : Numpy Array
        Factor used in nebular ionisation / dilute excitation calculations
        to account for the dilution of the radiation field.
    """
    outputs = ('w',)
    latex_name = ('W',)

class AtomicData(Input):
    outputs = ('atomic_data',)

class Abundance(Input):
    outputs = ('abundance',)

class Density(ArrayInput):
    outputs = ('density',)
    latex_name = ('\\rho',)

class TimeExplosion(Input):
    outputs = ('time_explosion',)
    latex_name = ('t_{\\textrm{exp}}',)

class JBlues(DataFrameInput):
    """
    Outputs:
    j_blues : Pandas DataFrame
        Mean intensity in the blue wing of each line.
    """
    outputs = ('j_blues',)
    latex_name = ('J_{lu}^{b}',)

class LinkTRadTElectron(Input):
    outputs = ('link_t_rad_t_electron',)
    latex_name = ('T_{\\textrm{electron}}/T_{\\textrm{rad}}',)
