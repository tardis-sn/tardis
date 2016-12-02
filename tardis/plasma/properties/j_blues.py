import numpy as np
import pandas as pd

from tardis.plasma.properties.base import ProcessingPlasmaProperty
from tardis.util import intensity_black_body


class JBluesBlackBody(ProcessingPlasmaProperty):
    '''
    Attributes
    ----------
    lte_j_blues : Pandas DataFrame, dtype float
                  J_blue values as calculated in LTE.
    '''
    outputs = ('j_blues',)
    latex_name = ('J^{b}_{lu(LTE)}')

    @staticmethod
    def calculate(lines, nu, t_rad):
        j_blues = intensity_black_body(nu.values[np.newaxis].T, t_rad)
        j_blues = pd.DataFrame(j_blues, index=lines.index,
                               columns=np.arange(len(t_rad)))
        return np.array(j_blues, copy=False)


class JBluesDiluteBlackBody(ProcessingPlasmaProperty):
    outputs = ('j_blues',)
    latex_name = ('J_{\\textrm{blue}}')

    @staticmethod
    def calculate(lines, nu, t_rad, w):
        j_blues = w * intensity_black_body(nu.values[np.newaxis].T, t_rad)
        j_blues = pd.DataFrame(j_blues, index=lines.index,
                               columns=np.arange(len(t_rad)))
        return np.array(j_blues, copy=False)
