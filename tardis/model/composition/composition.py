import pandas as pd
from astropy import units as u
from radioactivedecay import Inventory, Nuclide
from radioactivedecay.utils import Z_to_elem


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


class IsotopeMassFraction(pd.DataFrame):

    _metadata = ["time_0"]

    def __init__(self, *args, **kwargs):
        if "time_0" in kwargs:
            time_0 = kwargs["time_0"]
            kwargs.pop("time_0")
        else:
            time_0 = 0 * u.d
        super(IsotopeMassFraction, self).__init__(*args, **kwargs)
        self.time_0 = time_0

    @property
    def _constructor(self):
        return IsotopeMassFraction

    def _update_inventory(self):
        self.comp_dicts = [dict() for i in range(len(self.columns))]
        for (atomic_number, mass_number), abundances in self.iterrows():
            nuclear_symbol = f"{Z_to_elem(atomic_number)}{mass_number}"
            for i in range(len(self.columns)):
                self.comp_dicts[i][nuclear_symbol] = abundances[i]

    @classmethod
    def from_inventories(cls, inventories):
        multi_index_tuples = set([])
        for inventory in inventories:
            multi_index_tuples.update(
                [cls.id_to_tuple(key) for key in inventory.contents.keys()]
            )

        index = pd.MultiIndex.from_tuples(
            multi_index_tuples, names=["atomic_number", "mass_number"]
        )

        abundances = pd.DataFrame(
            data=0.0, index=index, columns=range(len(inventories))
        )

        for i, inventory in enumerate(inventories):
            for nuclide, abundance in inventory.masses("g").items():
                abundances.loc[cls.id_to_tuple(nuclide), i] = abundance

        return cls(abundances)

    @staticmethod
    def id_to_tuple(atomic_id):
        nuclide = Nuclide(atomic_id)
        return nuclide.Z, nuclide.A

    def to_inventories(self):
        """
        Convert DataFrame to a list of inventories interpreting the MultiIndex as
        atomic_number and mass_number

        Returns
        -------
        list
            list of radioactivedecay Inventories
        """

        comp_dicts = [dict() for i in range(len(self.columns))]
        for (atomic_number, mass_number), abundances in self.iterrows():
            nuclear_symbol = f"{Z_to_elem(atomic_number)}{mass_number}"
            for i in range(len(self.columns)):
                comp_dicts[i][nuclear_symbol] = abundances[i]
        return [Inventory(comp_dict, "g") for comp_dict in comp_dicts]

    def decay(self, t):
        """
        Decay the Model

        Parameters
        ----------
        t : float or astropy.units.Quantity
            if float it will be understood as days

        Returns
        -------
        pandas.DataFrame
            Decayed abundances
        """

        inventories = self.to_inventories()
        t_second = (
            u.Quantity(t, u.day).to(u.s).value - self.time_0.to(u.s).value
        )
        decayed_inventories = [item.decay(t_second) for item in inventories]
        df = IsotopeMassFraction.from_inventories(decayed_inventories)
        df.sort_index(inplace=True)
        return df

    def as_atoms(self):
        """
        Merge Isotope dataframe according to atomic number

        Returns
        -------
        pandas.DataFrame
            Merged isotope abundances
        """

        return self.groupby("atomic_number").sum()

    def merge(self, other, normalize=True):
        """
        Merge Isotope dataframe with abundance passed as parameter

        Parameters
        ----------
        other : pd.DataFrame
        normalize : bool
            If true, resultant dataframe will be normalized

        Returns
        -------
        pandas.DataFrame
            merged abundances
        """
        isotope_abundance = self.as_atoms()
        isotope_abundance = isotope_abundance.fillna(0.0)
        # Merge abundance and isotope dataframe
        modified_df = isotope_abundance.add(other, fill_value=0)

        if normalize:
            norm_factor = modified_df.sum(axis=0)
            modified_df /= norm_factor

        return modified_df
