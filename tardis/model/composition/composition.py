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
