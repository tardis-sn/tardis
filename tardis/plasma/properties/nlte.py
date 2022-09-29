import logging
import os

import numpy as np
import pandas as pd

from tardis.plasma.properties.base import (
    PreviousIterationProperty,
)
from tardis.plasma.properties.ion_population import PhiSahaNebular

__all__ = [
    "PreviousElectronDensities",
    "PreviousBetaSobolev",
]

logger = logging.getLogger(__name__)


class PreviousElectronDensities(PreviousIterationProperty):
    """
    Attributes
    ----------
    previous_electron_densities : The values for the electron densities converged upon in the previous iteration.
    """

    outputs = ("previous_electron_densities",)

    def set_initial_value(self, kwargs):
        initial_value = pd.Series(
            1000000.0,
            index=kwargs["abundance"].columns,
        )
        self._set_initial_value(initial_value)


class PreviousBetaSobolev(PreviousIterationProperty):
    """
    Attributes
    ----------
    previous_beta_sobolev : The beta sobolev values converged upon in the previous iteration.
    """

    outputs = ("previous_beta_sobolev",)

    def set_initial_value(self, kwargs):
        initial_value = pd.DataFrame(
            1.0,
            index=kwargs["atomic_data"].lines.index,
            columns=kwargs["abundance"].columns,
        )
        self._set_initial_value(initial_value)
