from tardis.iip_plasma.properties.base import ArrayInput, DataFrameInput, Input

__all__ = [
    "Abundance",
    "AtomicData",
    "Density",
    "DilutionFactor",
    "JBlues",
    "LinkTRadTElectron",
    "TRadiative",
    "TimeExplosion",
]


class TRadiative(ArrayInput):
    """
    Attributes
    ----------
    t_rad : Numpy Array, dtype float
    """

    outputs = ("t_rad",)
    latex_name = ("T_{\\textrm{rad}}",)


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


class Density(ArrayInput):
    """
    Attributes
    ----------
    density : Numpy array, dtype float
              Total density values
    """

    outputs = ("density",)
    latex_name = ("\\rho",)


class TimeExplosion(Input):
    """
    Attributes
    ----------
    ----------s
    time_explosion : Float
                     Time since explosion in seconds
    """

    outputs = ("time_explosion",)
    latex_name = ("t_{\\textrm{exp}}",)


class JBlues(DataFrameInput):
    """
    Attributes
    ----------
    j_blues : Pandas DataFrame
              Mean intensity in the blue wing of each line.
    """

    outputs = ("j_blues",)
    latex_name = ("J_{lu}^{b}",)


class LinkTRadTElectron(Input):
    """
    Attributes
    ----------
    link_t_rad_t_electron : Float
                            Value used for estimate of electron temperature.
                            Default is 0.9.
    """

    outputs = ("link_t_rad_t_electron",)
    latex_name = ("T_{\\textrm{electron}}/T_{\\textrm{rad}}",)
