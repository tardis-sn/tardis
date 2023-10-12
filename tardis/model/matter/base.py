from astropy import units as u


class Composition:
    """
    Holds information about model composition

    Parameters
    ----------
    density : astropy.units.quantity.Quantity
        An array of densities for each shell.
    isotopic_mass_fraction : pd.DataFrame
    atomic_mass : pd.DataFrame
    atomic_mass_unit: astropy.units.Unit

    Attributes
    ----------
    atomic_mass : pd.DataFrame
        Atomic mass of elements calculated for each shell.
    elemental_number_density : pd.DataFrame
        Number density of each element in each shell.
    """

    def __init__(
        self,
        density,
        elemental_mass_fraction,
        atomic_mass,
        atomic_mass_unit=u.g,
    ):
        self.density = density
        self.elemental_mass_fraction = elemental_mass_fraction
        self.atomic_mass_unit = atomic_mass_unit
        self._atomic_mass = atomic_mass

    @property
    def atomic_mass(self):
        """Atomic mass of elements in each shell"""
        if self._atomic_mass is None:
            raise AttributeError(
                "ModelState was not provided elemental masses."
            )
        return self._atomic_mass

    @property
    def elemental_number_density(self):
        """Elemental Number Density computed using the formula: (elemental_mass_fraction * density) / atomic mass"""
        if self.atomic_mass is None:
            raise AttributeError(
                "ModelState was not provided elemental masses."
            )
        return (self.elemental_mass_fraction * self.density).divide(
            self.atomic_mass, axis=0
        )


class ModelState:
    """
    Holds information about model geometry for radial 1D models.

    Parameters
    ----------
    composition : tardis.model.Composition
    geometry : tardis.model.geometry.radial1d.Radial1DGeometry
    time_explosion : astropy.units.quantity.Quantity

    Attributes
    ----------
    mass : pd.DataFrame
    number : pd.DataFrame
    """

    def __init__(self, composition, geometry, time_explosion):
        self.time_explosion = time_explosion
        self.composition = composition
        self.geometry = geometry

    @property
    def mass(self):
        """Mass calculated using the formula:
        mass_fraction * density * volume"""

        total_mass = (self.geometry.volume * self.composition.density).to(u.g)
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
