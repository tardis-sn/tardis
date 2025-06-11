import logging
from typing import Dict, Optional, Tuple

import pandas as pd
from astropy import units as u
from radioactivedecay import Inventory, Nuclide
from radioactivedecay.utils import Z_to_elem

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

    def to_inventories(self, cell_masses=None, cell_masses_unit="g"):
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
                if cell_masses is None:
                    comp_dicts[i][nuclear_symbol] = abundances[i]
                else:
                    comp_dicts[i][nuclear_symbol] = (
                        abundances[i] * cell_masses[i].to(u.g).value
                    )
        return [Inventory(comp_dict, cell_masses_unit) for comp_dict in comp_dicts]

    def calculate_decayed_mass_fractions(
        self, time_decay: u.Quantity, time_decay_unit: u.Unit = u.day
    ) -> "IsotopicMassFraction":
        """
        Decay the Model

        Parameters
        ----------
        time_decay : astropy.units.Quantity
            Time elapsed for decay calculation
        time_decay_unit : astropy.units.Unit, optional
            Unit for time_decay if it's provided as a scalar, by default u.day

        Returns
        -------
        IsotopicMassFraction
            Decayed abundances
        """
        inventories = self.to_inventories()
        time_decay = u.Quantity(time_decay, time_decay_unit)
        time_decay_second = time_decay.to(u.s).value - self.time_0.to(u.s).value
        logger.info("Decaying abundances for %.2g", time_decay)
        if time_decay_second < 0:
            logger.warning(
                "Decay time %f is negative. This could indicate a miss-specified input model."
                " A negative decay time can potentially lead to negative abundances.",
                time_decay_second,
            )
        decayed_inventories = [
            item.decay(time_decay_second, "s") for item in inventories
        ]
        df = IsotopicMassFraction.from_inventories(decayed_inventories)
        df = df.sort_index()
        assert (
            df.ge(0.0).all().all()
        ), "Negative abundances detected. Please make sure your input abundances are correct."
        return df

    def calculate_number_of_decays(
        self, time_decay: u.Quantity, cell_masses: Optional[u.Quantity] = None
    ) -> pd.DataFrame:
        """
        Calculate the number of decays over a given time period for each shell.

        Parameters
        ----------
        time_decay : ~astropy.units.Quantity
            Time elapsed over which to compute total decays. If provided without
            units, it is assumed to be in days.
        cell_masses : ~astropy.units.Quantity, optional
            1D array of shell masses. If provided, isotopic abundances will be
            multiplied by the corresponding shell mass (in grams).

        Returns
        -------
        pd.DataFrame
            Decays per isotope and per shell. Indexed by a MultiIndex
            (atomic_number, mass_number) with columns corresponding to shell indices.
        """
        # Convert time to days
        time_day = u.Quantity(time_decay, u.day).to(u.day).value

        # Build inventories for each shell (in grams if cell_masses given)
        inventories = self.to_inventories(
            cell_masses.to(u.g) if cell_masses is not None else None,
            cell_masses_unit="g",
        )

        decays_dict: Dict[Tuple[int, int, int], float] = {}
        for shell_idx, inventory in enumerate(inventories):
            shell_decays = inventory.cumulative_decays(time_day, "d")
            for isotope_sym, dec_count in shell_decays.items():
                nucl = Nuclide(isotope_sym)
                key = (nucl.Z, nucl.A, shell_idx)
                decays_dict[key] = dec_count

        # Build DataFrame
        decays_df = pd.DataFrame.from_dict(
            decays_dict, orient="index", columns=["decays"]
        )
        decays_df.index = pd.MultiIndex.from_tuples(
            decays_df.index, names=["atomic_number", "mass_number", "cell_id"]
        )
        decays_df = decays_df["decays"].unstack(level="cell_id").fillna(0.0)

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
