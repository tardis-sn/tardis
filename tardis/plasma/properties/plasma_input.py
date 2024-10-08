from tardis.plasma.properties.base import (
    Input,
    ProcessingPlasmaProperty,
)

__all__ = [
    "TRadiative",
    "DilutionFactor",
    "AtomicData",
    "Abundance",
    "NumberDensity",
    "IsotopeAbundance",
    "TimeExplosion",
    "JBlues",
    "LinkTRadTElectron",
    "HeliumTreatment",
    "ContinuumInteractionSpecies",
    "NLTEIonizationSpecies",
    "NLTEExcitationSpecies",
    "DilutePlanckianRadField",
]


class TRadiative(ProcessingPlasmaProperty):
    """
    Radiative temperature property.

    Attributes
    ----------
    t_rad : Numpy Array, dtype float
    """

    outputs = ("t_rad",)
    latex_name = (r"T_{\textrm{rad}}",)

    def calculate(self, dilute_planckian_radiation_field):
        return dilute_planckian_radiation_field.temperature.cgs.value


class DilutionFactor(ProcessingPlasmaProperty):
    """
    Dilution factor of the radiation field.

    Attributes
    ----------
    w : Numpy Array, dtype float between 0 and 1
        Factor used in nebular ionisation / dilute excitation calculations
        to account for the dilution of the radiation field.
    """

    outputs = ("w",)
    latex_name = ("W",)

    def calculate(self, dilute_planckian_radiation_field):
        return dilute_planckian_radiation_field.dilution_factor


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


class TimeExplosion(Input):
    """
    Attributes
    ----------
    time_explosion : Float
         Time since explosion in seconds
    """

    outputs = ("time_explosion",)
    latex_name = (r"t_{\textrm{exp}}",)


class JBlues(Input):
    """
    Attributes
    ----------
    j_blue_estimators : Numpy array
    """

    outputs = ("j_blues",)
    latex_name = (r"J_{\textrm{blue}}",)


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


class NLTEExcitationSpecies(Input):

    outputs = ("nlte_excitation_species",)


class NumberDensity(Input):
    """
    Attributes
    ----------
    number_density : Pandas DataFrame, dtype float
                     Indexed by atomic number, columns corresponding to zones
    """

    outputs = ("number_density",)
    latex_name = ("N_{i}",)


class DilutePlanckianRadField(Input):
    outputs = ("dilute_planckian_radiation_field",)
