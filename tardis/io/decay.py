import pandas as pd
from pyne import nucname, material
from astropy import units as u
from numpy import genfromtxt

class IsotopeAbundances(pd.DataFrame):

    @property
    def _constructor(self):
        return IsotopeAbundances

    def _update_material(self):
        self.comp_dicts = [{}] * len(self.columns)
        for (atomic_number, mass_number), abundances in self.iterrows():
            nuclear_symbol = '%s%d'.format(nucname.name(atomic_number),
                                           mass_number)
            for i in xrange(len(self.columns)):
                self.comp_dicts[i][nuclear_symbol] = abundances[i]

    @classmethod
    def from_materials(cls, materials):
        multi_index_tuples = set([])
        for material in materials:
            multi_index_tuples.update([cls.id_to_tuple(key)
                                       for key in material.keys()])

        index = pd.MultiIndex.from_tuples(
            multi_index_tuples, names=['atomic_number', 'mass_number'])


        abundances = pd.DataFrame(index=index, columns=xrange(len(materials)))

        for i, material in enumerate(materials):
            for key, value in material.items():
                abundances.loc[cls.id_to_tuple(key), i] = value

        return abundances




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

        comp_dicts = [{}] * len(self.columns)
        for (atomic_number, mass_number), abundances in self.iterrows():
            nuclear_symbol = '{0:s}{1:d}'.format(nucname.name(atomic_number),
                                           mass_number)
            for i in xrange(len(self.columns)):
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
        t_second = u.Quantity(t, u.day).to(u.s).value

        decayed_materials = [item.decay(t_second) for item in materials]

        df = IsotopeAbundances.from_materials(decayed_materials)
        df.sort_index(inplace=True)
        return df

def read_isotope_file(fname):
    data = genfromtxt(fname, dtype=None, skip_header=1)
    ids = [nucname.id(i) for i, _ in data]
    atomic_number = [nucname.znum(i) for i in ids]
    mass_number = [nucname.anum(i) for i in ids]

    mass = [i for _, i in data]
    combined = zip(atomic_number, mass_number)
    index = pd.MultiIndex.from_tuples(combined,
                                      names=['atomic_number', 'mass_number'])
    df = pd.DataFrame(mass, index=index, columns=['mass'])
    return df