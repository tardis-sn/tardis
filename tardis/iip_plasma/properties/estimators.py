from tardis.iip_plasma.properties.base import Input

__all__ = [
    "BfHeatingEstimator",
    "CollDeexcHeatingEstimator",
    "FfHeatingEstimator",
    "PhotoIonRateEstimator",
    "PhotoIonRateStatistics",
    "StimRecombCoolingEstimator",
    "StimRecombRateEstimator",
]


class PhotoIonRateEstimator(Input):
    """
    Attributes
    ----------
    photo_ion_estimator : Pandas DataFrame
                          Monte Carlo estimator for the photoionization rate coefficient.
    """

    outputs = ("photo_ion_estimator",)
    latex_name = ("\\gamma_{\\textrm{estimator}}",)


class PhotoIonRateStatistics(Input):
    """
    Attributes
    ----------
    photo_ion_statistics : Pandas DataFrame
                           Number of updates of the photoionization rate estimator.
    """

    outputs = ("photo_ion_statistics",)
    latex_name = ("\\N_{\\textrm{bf_estimator}}",)


class StimRecombRateEstimator(Input):
    """
    Attributes
    ----------
    stim_recomb_estimator : Pandas DataFrame
                            Monte Carlo estimator for the stimulated recombination rate coefficient.
    """

    outputs = ("stim_recomb_estimator",)
    latex_name = ("\\alpha_{\\textrm{estimator}}^{\\textrm{st}}",)


class BfHeatingEstimator(Input):
    """
    Attributes
    ----------
    bf_heating_estimator : Pandas DataFrame
                           Monte Carlo estimator for the bound-free heating rate coefficient.
    """

    outputs = ("bf_heating_estimator",)
    latex_name = ("\\gamma_{\\textrm{estimator}}^{\\textrm{E}}",)


class StimRecombCoolingEstimator(Input):
    """
    Attributes
    ----------
    stim_recomb_cooling_estimator : Pandas DataFrame
                           Monte Carlo estimator for the stimulated recombination cooling rate coefficient.
    """

    outputs = ("stim_recomb_cooling_estimator",)
    latex_name = ("\\alpha_{\\textrm{estimator}}^{\\textrm{E,st}}",)


class FfHeatingEstimator(Input):
    """
    Attributes
    ----------
    ff_heating_estimator : Pandas DataFrame
                           Monte Carlo estimator for the free-free heating rate coefficient.
    """

    outputs = ("ff_heating_estimator",)
    latex_name = ("\\h_{\\textrm{estimator}}^{\\textrm{ff}}",)


class CollDeexcHeatingEstimator(Input):
    """
    Attributes
    ----------
    coll_deexc_heating_estimator : Pandas DataFrame
                           Monte Carlo estimator for the collisional deexcitation heating rate coefficient.
    """

    outputs = ("coll_deexc_heating_estimator",)
    latex_name = ("\\h_{\\textrm{estimator}}^{\\textrm{coll_deexc}}",)
