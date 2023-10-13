from astropy import units as u
from tardis.model.matter.decay import IsotopeAbundances


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
        nuclide_mass_fraction,
        atomic_mass,
        atomic_mass_unit=u.g,
    ):
        self.density = density
        self.nuclide_mass_fraction = nuclide_mass_fraction
        # self.elemental_mass_fraction = elemental_mass_fraction
        self.atomic_mass_unit = atomic_mass_unit
        self._atomic_mass = atomic_mass
        """
        self.raw_abundance = self._abundance
        self.raw_isotope_abundance = isotope_abundance
        
        if elemental_mass is not None:
            mass = {}
            stable_atomic_numbers = self.raw_abundance.index.to_list()
            for z in stable_atomic_numbers:
                mass[z] = [
                    elemental_mass[z]
                    for i in range(self.raw_abundance.columns.size)
                ]
            stable_isotope_mass = pd.DataFrame(mass).T

            isotope_mass = {}
            for atomic_number, i in self.raw_isotope_abundance.decay(
                self.time_explosion
            ).groupby(level=0):
                i = i.loc[atomic_number]
                for column in i:
                    mass = {}
                    shell_abundances = i[column]
                    isotopic_masses = [
                        rd.Nuclide(Z_DICT[atomic_number] + str(i)).atomic_mass
                        for i in shell_abundances.index.to_numpy()
                    ]
                    mass[atomic_number] = (
                        shell_abundances * isotopic_masses
                    ).sum()
                    mass[atomic_number] /= shell_abundances.sum()
                    mass[atomic_number] = mass[atomic_number] * u.u.to(u.g)
                    if isotope_mass.get(column) is None:
                        isotope_mass[column] = {}
                    isotope_mass[column][atomic_number] = mass[atomic_number]
            isotope_mass = pd.DataFrame(isotope_mass)

            atomic_mass = pd.concat([stable_isotope_mass, isotope_mass])
        """

    @property
    def isotope_mass_fraction(self):
        filtered_nuclide_mass_fraction = self.nuclide_mass_fraction[
            self.nuclide_mass_fraction.index.get_level_values(1) != -1
        ]
        return IsotopeAbundances(filtered_nuclide_mass_fraction)

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


class MatterState:
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
