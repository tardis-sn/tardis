import numpy as np
import pandas as pd
import radioactivedecay as rd
from astropy import units as u
from radioactivedecay.decaydata import DEFAULTDATA as RD_DEFAULT_DATA

from tardis.model.matter.decay import IsotopicMassFraction


def compile_rd_isotope_masses():
    """
    Compiles the masses of isotopes from the default data in RD_DEFAULT_DATA.

    Parameters
    ----------
    None

    Returns
    -------
    pandas.Series
        A series containing the masses of isotopes, indexed by atomic number and mass number.
    """
    atomic_numbers = []
    mass_numbers = []
    nuclide_masses = []
    for nuclide_name in RD_DEFAULT_DATA.nuclides:
        current_nuclide = rd.Nuclide(nuclide_name)
        atomic_numbers.append(current_nuclide.Z)
        mass_numbers.append(current_nuclide.A)
        nuclide_masses.append(current_nuclide.atomic_mass)

    isotope_mass_index = pd.MultiIndex.from_arrays(
        (atomic_numbers, mass_numbers), names=["atomic_number", "mass_number"]
    )
    isotope_masses = pd.Series(
        index=isotope_mass_index, data=nuclide_masses
    ).sort_index()
    # there are duplicates that are likely due to excited states
    # dropping them for now

    return isotope_masses.drop_duplicates()


ISOTOPE_MASSES = compile_rd_isotope_masses()


class Composition:
    """
    Holds information about model composition

    Parameters
    ----------
    density : astropy.units.quantity.Quantity
        An array of densities for each shell.
    isotopic_mass_fraction : pd.DataFrame
    raw_isotope_abundance : pd.DataFrame
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
        raw_isotope_abundance,
        element_masses,
        element_masses_unit=u.g,
    ):
        self.density = density
        assert np.all(
            nuclide_mass_fraction.values >= 0
        ), "Negative mass fraction detected"
        self.nuclide_mass_fraction = nuclide_mass_fraction

        self.nuclide_masses_unit = element_masses_unit

        self.nuclide_masses = element_masses
        self.nuclide_masses.index = self.convert_element2nuclide_index(
            element_masses.index
        )

        isotope_masses = self.assemble_isotope_masses()

        self.nuclide_masses = pd.concat([self.nuclide_masses, isotope_masses])
        self.raw_isotope_abundance = raw_isotope_abundance

    def assemble_isotope_masses(self):
        isotope_mass_df = pd.Series(
            index=self.isotopic_mass_fraction.index, data=-1
        )
        for isotope_tuple in self.isotopic_mass_fraction.index:
            isotope_symbol = int("{:03d}{:03d}0000".format(*isotope_tuple))
            isotope_mass = rd.Nuclide(isotope_symbol).atomic_mass * u.u.to(u.g)
            isotope_mass_df[isotope_tuple] = isotope_mass

        return isotope_mass_df

    @staticmethod
    def convert_element2nuclide_index(element_index):
        new_nuclide_index = pd.MultiIndex.from_product([element_index, [-1]])
        new_nuclide_index.names = ["atomic_number", "mass_number"]
        return new_nuclide_index

    @property
    def isotopic_mass_fraction(self):
        filtered_nuclide_mass_fraction = self.nuclide_mass_fraction[
            self.nuclide_mass_fraction.index.get_level_values(1) != -1
        ]
        return IsotopicMassFraction(filtered_nuclide_mass_fraction)

    @property
    def elemental_mass_fraction(self):
        return self.nuclide_mass_fraction.groupby(level=0).sum()

    @property
    def element_masses(self):
        """Atomic mass of elements in each shell"""
        element_masses = self.nuclide_masses[
            self.nuclide_masses.index.get_level_values(1) == -1
        ]
        element_masses.index = element_masses.index.droplevel(1)
        return element_masses

    @property
    def effective_element_masses(self):
        # This is functionality that we will likely want to remove
        effective_element_masses = self.nuclide_mass_fraction[
            self.nuclide_mass_fraction.index.get_level_values(1) == -1
        ].copy()
        effective_element_masses.index = (
            effective_element_masses.index.droplevel(1)
        )
        for col in effective_element_masses.columns:
            effective_element_masses[col] = self.element_masses.loc[
                effective_element_masses.index
            ]

        current_isotope_masses = ISOTOPE_MASSES.loc[
            self.isotopic_mass_fraction.index
        ]
        contributing_isotope_masses = (
            self.isotopic_mass_fraction.multiply(current_isotope_masses, axis=0)
            .groupby(level=0)
            .sum()
        )
        effective_isotopes_masses = (
            contributing_isotope_masses
            / self.isotopic_mass_fraction.groupby(level=0).sum()
        ) * u.u.to(u.g)

        effective_element_masses = pd.concat(
            [effective_element_masses, effective_isotopes_masses]
        )
        return effective_element_masses

    @property
    def elemental_number_density(self):
        """Elemental Number Density computed using the formula: (elemental_mass_fraction * density) / atomic mass"""
        return (
            self.elemental_mass_fraction
            * self.density.to(u.g / u.cm**3).value
        ).divide(
            self.effective_element_masses.reindex(
                self.elemental_mass_fraction.index
            ),
            axis=0,
        )

    @property
    def isotopic_number_density(self):
        """Isotopic Number Density computed using the formula: (isotopic_mass_fraction * density) / atomic mass"""
        return (
            self.isotopic_mass_fraction * self.density.to(u.g / u.cm**3).value
        ).divide(
            ISOTOPE_MASSES.loc[self.isotopic_mass_fraction.index] * u.u.to(u.g),
            axis=0,
        )

    def calculate_mass_fraction_at_time(self, time_explosion):
        """
        Calculate the mass fraction at a given time using the radioactive decay from
        the IsotopicMassFraction.

        Parameters
        ----------
        time_explosion : astropy.units.quantity.Quantity
            The time of the explosion.

        Returns
        -------
        None

        Examples
        --------
        >>> composition.calculate_mass_fraction_at_time(10 * u.s)
        """
        if self.isotopic_mass_fraction.empty:
            return self.elemental_mass_fraction
        else:
            self.isotopic_mass_fraction.decay(time_explosion)

    def calculate_elemental_cell_masses(self, volume):
        """
        Calculate the elemental cell masses.

        Parameters
        ----------
        volume : astropy.units.quantity.Quantity
        The volume of the cell.

        Returns
        -------
        numpy.ndarray
        An array of elemental cell masses.

        Examples
        --------
        >>> composition.calculate_cell_masses(10 * u.cm**3)
        """
        return (
            self.elemental_mass_fraction * (self.density * volume).to(u.g).value
        )

    def calculate_cell_masses(self, volume):
        """
        Calculate the cell masses.

        Parameters
        ----------
        volume : astropy.units.quantity.Quantity
        The volume of the cell.

        Returns
        -------
        astropy.units.quantity.Quantity
        An array of cell masses.

        Examples
        --------
        >>> composition.calculate_cell_masses(10 * u.cm**3)
        """
        return (self.density * volume).to(u.g)
