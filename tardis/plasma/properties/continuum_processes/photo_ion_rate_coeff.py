from tardis.plasma.properties.base import Input
from tardis.plasma.properties.continuum_processes.rates import H


class PhotoIonRateCoeff(Input):
    """
    Attributes
    ----------
    gamma_estimator : pandas.DataFrame, dtype float
        Unnormalized MC estimator for the rate coefficient for radiative
        ionization.
    """

    outputs = ("gamma",)
    latex_name = (r"\gamma",)
