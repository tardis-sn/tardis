import pandas as pd

from tardis.plasma.properties.base import (
    Input,
    ProcessingPlasmaProperty,
)

__all__ = [
    "Abundance",
    "AtomicData",
    "ContinuumInteractionSpecies",
    "DilutePlanckianRadField",
    "DilutionFactor",
    "ElectronDensitiesInput",
    "HeliumTreatment",
    "IonNumberDensityInput",
    "IsotopeAbundance",
    "JBlues",
    "LevelNumberDensityInput",
    "LinkTRadTElectron",
    "NumberDensity",
    "TRadiative",
    "TimeExplosion",
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

    def calculate(self, dilute_planckian_radiation_field: object) -> object:
        """Return the radiation-field temperature in cgs units."""
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

    def calculate(self, dilute_planckian_radiation_field: object) -> object:
        """Return the radiation-field dilution factor."""
        return dilute_planckian_radiation_field.dilution_factor


class AtomicData(Input):
    """Store the atomic-data input."""

    outputs = ("atomic_data",)


class Abundance(Input):
    """Store elemental abundances."""

    outputs = ("abundance",)


class IsotopeAbundance(Input):
    """Store isotope abundances."""

    outputs = ("isotope_abundance",)


class TimeExplosion(Input):
    """Store time since explosion in seconds."""

    outputs = ("time_explosion",)
    latex_name = (r"t_{\textrm{exp}}",)


class JBlues(Input):
    """Store line blue-wing mean intensities."""

    outputs = ("j_blues",)
    latex_name = (r"J_{\textrm{blue}}",)


class ElectronDensitiesInput(Input):
    """Electron densities accepted from an external equilibrium solve."""

    electron_densities: pd.Series
    outputs = ("electron_densities",)


class IonNumberDensityInput(Input):
    """Ion number densities accepted from an external equilibrium solve."""

    ion_number_density: pd.DataFrame
    outputs = ("ion_number_density",)


class LevelNumberDensityInput(Input):
    """Level number densities accepted from an external equilibrium solve."""

    level_number_density: pd.DataFrame
    outputs = ("level_number_density",)


class LinkTRadTElectron(Input):
    """Store the electron-to-radiation temperature ratio."""

    outputs = ("link_t_rad_t_electron",)
    latex_name = (r"T_{\textrm{electron}}/T_{\textrm{rad}}",)


class HeliumTreatment(Input):
    """Store the configured helium treatment."""

    outputs = ("helium_treatment",)


class ContinuumInteractionSpecies(Input):
    """Store species with enabled continuum interactions."""

    outputs = ("continuum_interaction_species",)


class NumberDensity(Input):
    """Store elemental number densities by shell."""

    outputs = ("number_density",)
    latex_name = ("N_{i}",)


class DilutePlanckianRadField(Input):
    """Store the dilute Planckian radiation field."""

    outputs = ("dilute_planckian_radiation_field",)
