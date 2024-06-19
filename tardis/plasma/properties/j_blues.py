import numpy as np
import pandas as pd

from tardis import constants as const
from tardis.plasma.properties.base import (
    DataFrameInput,
    ProcessingPlasmaProperty,
)
from tardis.util.base import intensity_black_body


class JBluesBlackBody(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    lte_j_blues : Pandas DataFrame, dtype float
                  J_blue values as calculated in LTE.
    """

    outputs = ("j_blues",)
    latex_name = "J^{b}_{lu(LTE)}"

    @staticmethod
    def calculate(lines, nu, t_rad):
        j_blues = intensity_black_body(nu.values[np.newaxis].T, t_rad)
        j_blues = pd.DataFrame(
            j_blues, index=lines.index, columns=np.arange(len(t_rad))
        )
        return j_blues


class JBluesDiluteBlackBody(ProcessingPlasmaProperty):
    outputs = ("j_blues",)
    latex_name = r"J_{\textrm{blue}}"

    @staticmethod
    def calculate(lines, nu, t_rad, w):
        j_blues = w * intensity_black_body(nu.values[np.newaxis].T, t_rad)
        j_blues = pd.DataFrame(
            j_blues, index=lines.index, columns=np.arange(len(t_rad))
        )
        return j_blues


class JBluesDetailed(ProcessingPlasmaProperty):
    outputs = ("j_blues",)
    latex_name = "J_{\\textrm{blue}}"

    def __init__(self, plasma_parent, w_epsilon):
        super(JBluesDetailed, self).__init__(plasma_parent)
        self.w_epsilon = w_epsilon

    def calculate(
        self, lines, nu, t_rad, w, j_blues_norm_factor, j_blue_estimator
    ):
        # Used for initialization
        if len(j_blue_estimator) == 0:
            return JBluesDiluteBlackBody.calculate(lines, nu, t_rad, w)
        else:
            j_blues = pd.DataFrame(
                j_blue_estimator * j_blues_norm_factor.value,
                index=lines.index,
                columns=np.arange(len(t_rad)),
            )

            for i in range(len(t_rad)):
                zero_j_blues = j_blues[i] == 0.0
                j_blues[i][
                    zero_j_blues
                ] = self.w_epsilon * intensity_black_body(
                    nu[zero_j_blues].values, t_rad[i]
                )
            return j_blues


class JBluesNormFactor(ProcessingPlasmaProperty):

    outputs = ("j_blues_norm_factor",)
    latex = (
        r"\frac{c time_\textrm{simulation}}}{4\pi"
        r"time_\textrm{simulation} volume}"
    )

    @staticmethod
    def calculate(time_explosion, time_simulation, volume):
        return (
            const.c.cgs
            * time_explosion
            / (4 * np.pi * time_simulation * volume)
        )


class JBluesEstimator(DataFrameInput):
    outputs = ("j_blue_estimator",)
