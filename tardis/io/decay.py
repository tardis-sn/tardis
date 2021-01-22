import pandas as pd
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
            nuclear_symbol = f"{nucname.name(atomic_number)}{mass_number}"
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

    @staticmethod
    def id_to_tuple(atomic_id):
        return nucname.znum(atomic_id), nucname.anum(atomic_id)

    def to_materials(self):
        """
        Convert DataFrame to a list of materials interpreting the MultiIndex as
        atomic_number and mass_number

        Returns
        -------
            : ~list
            list of pyne Materialss
        :return:
        """

        comp_dicts = [dict() for i in range(len(self.columns))]
        for (atomic_number, mass_number), abundances in self.iterrows():
            nuclear_symbol = "{0:s}{1:d}".format(
                nucname.name(atomic_number), mass_number
            )
            for i in range(len(self.columns)):
                comp_dicts[i][nuclear_symbol] = abundances[i]
        return [material.Material(comp_dict) for comp_dict in comp_dicts]

    def decay(self, t):
        """
        Decay the Model

        Parameters
        ----------

        t: ~float or ~astropy.units.Quantity
            if float it will be understood as days

        Returns:
            : decayed abundances
        """

        materials = self.to_materials()
        t_second = (
            u.Quantity(t, u.day).to(u.s).value - self.time_0.to(u.s).value
        )
        decayed_materials = [item.decay(t_second) for item in materials]
        for i in range(len(materials)):
            materials[i].update(decayed_materials[i])
        df = IsotopeAbundances.from_materials(materials)
        df.sort_index(inplace=True)
        return df

    def as_atoms(self):
        """
        Merge Isotope dataframe according to atomic number

        Returns:
            : merged isotope abundances
        """

        return self.groupby("atomic_number").sum()

    def merge(self, other, normalize=True):
        """
        Merge Isotope dataframe with abundance passed as parameter

        Parameters
        ----------
        other: pd.DataFrame
        normalize : bool
            If true, resultant dataframe will be normalized

        Returns:
            : merged abundances
        """
        isotope_abundance = self.as_atoms()
        isotope_abundance = isotope_abundance.fillna(0.0)
        # Merge abundance and isotope dataframe
        modified_df = isotope_abundance.add(other, fill_value=0)

        if normalize:
            norm_factor = modified_df.sum(axis=0)
            modified_df /= norm_factor

        return modified_df
