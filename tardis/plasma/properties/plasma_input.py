from tardis.plasma.properties.base import Input, ArrayInput, DataFrameInput

__all__ = [
    "TRadiative",
    "DilutionFactor",
    "AtomicData",
    "Abundance",
    "IsotopeAbundance",
    "Density",
    "TimeExplosion",
    "JBlueEstimator",
    "LinkTRadTElectron",
    "HeliumTreatment",
    "RInner",
    "TInner",
    "Volume",
    "ContinuumInteractionSpecies",
    "NLTEIonizationSpecies",
]


class TRadiative(ArrayInput):
    """
    Attributes
    ----------
    t_rad : Numpy Array, dtype float
    """

    outputs = ("t_rad",)
    latex_name = (r"T_{\textrm{rad}}",)


class DilutionFactor(ArrayInput):
    """
    Attributes
    ----------
    w : Numpy Array, dtype float between 0 and 1
        Factor used in nebular ionisation / dilute excitation calculations
        to account for the dilution of the radiation field.
    """

    outputs = ("w",)
    latex_name = ("W",)


class AtomicData(Input):
    """
    Attributes
    ----------
    atomic_data : Object
    """

    outputs = ("atomic_data",)


class Abundance(Input):
    """
    Attributes
    ----------
    abundance : Numpy array, dtype float
        Fractional abundance of elements
    """

    outputs = ("abundance",)


class IsotopeAbundance(Input):
    """
    Attributes
    ----------
    isotope_abundance : Numpy array, dtype float
        Fractional abundance of isotopes
    """

    outputs = ("isotope_abundance",)


class Density(ArrayInput):
    """
    Attributes
    ----------
    density : Numpy array, dtype float
      Total density values
    """

    outputs = ("density",)
    latex_name = (r"\rho",)


class TimeExplosion(Input):
    """
    Attributes
    ----------
    time_explosion : Float
         Time since explosion in seconds
    """

    outputs = ("time_explosion",)
    latex_name = (r"t_{\textrm{exp}}",)


class JBlueEstimator(ArrayInput):
    """
    Attributes
    ----------
    j_blue_estimators : Numpy array
    """

    outputs = ("j_blue_estimators",)
    latex_name = (r"J_{\textrm{blue-estimator}}",)


class LinkTRadTElectron(Input):
    """
    Attributes
    ----------
    link_t_rad_t_electron : Float
        Value used for estimate of electron temperature.
        Default is 0.9.
    """

    outputs = ("link_t_rad_t_electron",)
    latex_name = (r"T_{\textrm{electron}}/T_{\textrm{rad}}",)


class HeliumTreatment(Input):
    outputs = ("helium_treatment",)


class RInner(Input):
    outputs = ("r_inner",)


class TInner(Input):
    outputs = ("t_inner",)


class Volume(Input):
    outputs = ("volume",)


class ContinuumInteractionSpecies(Input):
    """
    Attributes
    ----------
    continuum_interaction_species : Pandas MultiIndex, dtype int
        Atomic and ion numbers of elements for which continuum interactions
        (radiative/collisional ionization and recombination) are treated
    """

    outputs = ("continuum_interaction_species",)


class NLTEIonizationSpecies(Input):

    outputs = ("nlte_ionization_species",)
