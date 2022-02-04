import pandas as pd
from radioactivedecay import Nuclide, Inventory
from radioactivedecay.utils import Z_to_elem
from pyne import nucname, material
from astropy import units as u


class IsotopeAbundances(pd.DataFrame):

    _metadata = ["time_0"]

    def __init__(self, *args, **kwargs):
        if "time_0" in kwargs:
            time_0 = kwargs["time_0"]
            kwargs.pop("time_0")
        else:
            time_0 = 0 * u.d
        super(IsotopeAbundances, self).__init__(*args, **kwargs)
        self.time_0 = time_0

    @property
    def _constructor(self):
        return IsotopeAbundances

    def _update_material(self):
        self.comp_dicts = [dict() for i in range(len(self.columns))]
        for (atomic_number, mass_number), abundances in self.iterrows():
            nuclear_symbol = f"{Z_to_elem(atomic_number)}{mass_number}"
            for i in range(len(self.columns)):
                self.comp_dicts[i][nuclear_symbol] = abundances[i]

    @classmethod
    def from_materials(cls, materials):
        multi_index_tuples = set([])
        for material in materials:
            multi_index_tuples.update(
                [cls.id_to_tuple(key) for key in material.keys()]
            )

        index = pd.MultiIndex.from_tuples(
            multi_index_tuples, names=["atomic_number", "mass_number"]
        )

        abundances = pd.DataFrame(
            data=0.0, index=index, columns=range(len(materials))
        )

        for i, material in enumerate(materials):
            for key, value in material.items():
                abundances.loc[cls.id_to_tuple(key), i] = value

        return cls(abundances)
    
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
            for nuclide, abundance in inventory.mass_fractions():
                abundances.loc[cls.id_to_tuple(nuclide), i] = abundance

        return cls(abundances)

    @staticmethod
    def id_to_tuple(atomic_id):
        nuclide = Nuclide(atomic_id)
        return nuclide.Z, nuclide.A

    def to_materials(self):
        """
        Convert DataFrame to a list of materials interpreting the MultiIndex as
        atomic_number and mass_number

        Returns
        -------
        list
            list of pyne Materials
        """

        comp_dicts = [dict() for i in range(len(self.columns))]
        for (atomic_number, mass_number), abundances in self.iterrows():
            nuclear_symbol = f"{Z_to_elem(atomic_number)}{mass_number}"
            for i in range(len(self.columns)):
                comp_dicts[i][nuclear_symbol] = abundances[i]
        return [material.Material(comp_dict) for comp_dict in comp_dicts]
    
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
        return [Inventory(comp_dict, "kg") for comp_dict in comp_dicts]

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
        for i in range(len(inventories)):
            inventories[i].update(decayed_inventories[i])
        df = IsotopeAbundances.from_inventories(inventories)
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
