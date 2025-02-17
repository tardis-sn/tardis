import logging

import pandas as pd
from astropy import units as u
from radioactivedecay import Inventory, Nuclide
from radioactivedecay.utils import Z_to_elem, elem_to_Z

logger = logging.getLogger(__name__)


class IsotopicMassFraction(pd.DataFrame):
    _metadata = ["time_0"]

    def __init__(self, *args, **kwargs):
        if "time_0" in kwargs:
            time_0 = kwargs["time_0"]
            kwargs.pop("time_0")
        else:
            time_0 = 0 * u.d
        super().__init__(*args, **kwargs)
        self.time_0 = time_0

    @property
    def _constructor(self):
        return IsotopicMassFraction

    def _update_inventory(self):
        self.comp_dicts = [{} for i in range(len(self.columns))]
        for (atomic_number, mass_number), mass_fractions in self.iterrows():
            nuclear_symbol = f"{Z_to_elem(atomic_number)}{mass_number}"
            for i in range(len(self.columns)):
                self.comp_dicts[i][nuclear_symbol] = mass_fractions[i]

    @classmethod
    def from_inventories(cls, inventories):
        multi_index_tuples = set()
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

    def to_inventories(self, shell_masses=None):
        """
        Convert DataFrame to a list of inventories interpreting the MultiIndex as
        atomic_number and mass_number

        Returns
        -------
        list
            list of radioactivedecay Inventories
        """
        comp_dicts = [{} for i in range(len(self.columns))]
        for (atomic_number, mass_number), abundances in self.iterrows():
            nuclear_symbol = f"{Z_to_elem(atomic_number)}{mass_number}"
            for i in range(len(self.columns)):
                if shell_masses is None:
                    comp_dicts[i][nuclear_symbol] = abundances[i]
                else:
                    comp_dicts[i][nuclear_symbol] = (
                        abundances[i] * shell_masses[i].to(u.g).value
                    )
        return [Inventory(comp_dict, "g") for comp_dict in comp_dicts]

    def calculate_decayed_mass_fractions(self, time_decay):
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
            u.Quantity(time_decay, u.day).to(u.s).value - self.time_0.to(u.s).value
        )
        logger.info(f"Decaying abundances for {t_second} seconds")
        if t_second < 0:
            logger.warning(
                "Decay time %f is negative. This could indicate a miss-specified input model."
                " A negative decay time can potentially lead to negative abundances.",
                t_second,
            )
        decayed_inventories = [item.decay(t_second, "s") for item in inventories]
        df = IsotopicMassFraction.from_inventories(decayed_inventories)
        df = df.sort_index()
        assert (
            df.ge(0.0).all().all()
        ), "Negative abundances detected. Please make sure your input abundances are correct."
        return df

    def calculate_number_of_decays(self, time_decay, shell_masses=None):
        """
        Calculate the number of decays over a given time period for each shell.

        Parameters
        ----------
        time_decay : astropy.units.Quantity
            Time elapsed for which to compute the total decays (converted to days internally).
        shell_masses : astropy.units.Quantity, optional
            Masses of each shell. If provided, isotopes in each shell's inventory will be
            multiplied by the corresponding shell mass.

        Returns
        -------
        pandas.DataFrame
            A DataFrame of decays indexed by (atomic_number, mass_number),
            with columns representing the decays in each shell.
        """
        # Convert the time to days for radioactivedecay
        t_days = u.Quantity(time_decay, u.day).value

        # Get inventories for each shell. If shell_masses is provided,
        # each inventory is in grams of isotopes.
        inventories = self.to_inventories(shell_masses)

        decays_dict = {}
        for shell_index, inventory in enumerate(inventories):
            shell_decays = inventory.cumulative_decays(t_days, "d")
            for isotope, dec_count in shell_decays.items():
                element_symbol, massno_str = isotope.split("-")
                atomic_number = elem_to_Z(element_symbol)
                mass_number = int(massno_str)
                decays_dict[(atomic_number, mass_number, shell_index)] = dec_count

        # Build a decays DataFrame
        decays_df = pd.DataFrame.from_dict(
            decays_dict, orient="index", columns=["decays"]
        )
        decays_df.index = pd.MultiIndex.from_tuples(
            decays_df.index, names=["atomic_number", "mass_number", "shell"]
        )
        decays_df = decays_df["decays"].unstack(level="shell").fillna(0)

        return decays_df

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
