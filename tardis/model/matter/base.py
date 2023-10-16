from astropy import units as u


class MatterState:
    """
    Holds information about model geometry for radial 1D models.

    Parameters
    ----------
    composition : tardis.model.Composition

    Attributes
    ----------
    mass : pd.DataFrame
    number : pd.DataFrame
    """

    def __init__(self, composition, time_explosion):
        self.time_explosion = time_explosion
        self.composition = composition
        self.geometry = geometry

    def calculate_mass(self, geometry):
        """Mass calculated using the formula:
        mass_fraction * density * volume"""

        total_mass = (geometry.volume * self.composition.density).to(u.g)
        return self.composition.elemental_mass_fraction * total_mass.value

    @property
    def number(self):
        """Number calculated using the formula:
        mass / atomic_mass"""
        if self.composition.atomic_mass is None:
            raise AttributeError(
                "ModelState was not provided elemental masses."
            )
        return (self.mass).divide(self.composition.atomic_mass, axis=0)
