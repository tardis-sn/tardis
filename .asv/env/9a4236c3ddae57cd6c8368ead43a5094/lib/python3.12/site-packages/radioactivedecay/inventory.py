"""
The inventory module defines the ``Inventory`` and ``InventoryHP`` classes. Instances of these
classes contain one or more nuclides each with an associated amount (number of atoms). The decay
of the nuclide(s) can be calculated by using the ``decay()`` method. The ``Inventory`` class
performs calculations using double-precision floats. The ``InventoryHP`` class performs
calculations using SymPy high numerical precision operations. A ``DecayData`` dataset is associated
with ``Inventory`` and ``InventoryHP`` instances (default is rd.DEFAULTDATA).

The docstring code examples assume that ``radioactivedecay`` has been imported
as ``rd``:

.. highlight:: python
.. code-block:: python

    >>> import radioactivedecay as rd

"""

import collections
import csv
import itertools
import numbers
import pathlib
from abc import ABC, abstractmethod
from typing import Any, Callable, Dict, List, Optional, Tuple, Type, Union

import matplotlib
import numpy as np
import pandas as pd
from scipy import sparse
from sympy import Integer, Matrix, exp, nsimplify
from sympy.core.expr import Expr

from radioactivedecay.converters import (
    QuantityConverter,
    QuantityConverterFloat,
    QuantityConverterSympy,
    UnitConverter,
    UnitConverterFloat,
    UnitConverterSympy,
)
from radioactivedecay.decaydata import (
    DEFAULTDATA,
    DecayData,
    DecayMatrices,
    DecayMatricesScipy,
    DecayMatricesSympy,
)
from radioactivedecay.nuclide import Nuclide
from radioactivedecay.plots import decay_graph
from radioactivedecay.utils import (
    add_dictionaries,
    parse_nuclide,
    sort_dictionary_alphabetically,
    sort_list_according_to_dataset,
)

# functools.singledispatchmethod is new in Python version 3.8
# this block adds support manually for Python <3.8
# the except block can be removed once min required Python version is 3.8
try:
    from functools import singledispatchmethod
except ImportError:
    from functools import singledispatch, update_wrapper

    def singledispatchmethod(func: Callable[..., Any]):  # type:ignore
        """Adds singledispatch support for instance methods of a class."""

        dispatcher = singledispatch(func)

        def wrapper(*args, **kwargs):  # type:ignore
            return dispatcher.dispatch(args[1].__class__)(*args, **kwargs)

        wrapper.register = dispatcher.register
        update_wrapper(wrapper, func)
        return wrapper


def _write_csv_file(
    filename: Union[str, pathlib.Path],
    rows: List[List[str]],
    delimiter: str,
    encoding: str,
) -> None:
    """Write a CSV file from a list of rows."""

    with open(filename, "w", encoding=encoding, newline="") as file:
        writer_object = csv.writer(file, delimiter=delimiter)
        writer_object.writerows(rows)


# pylint: disable=too-many-arguments, too-many-lines, too-many-locals


class AbstractInventory(ABC):
    """
    Template class for float (SciPy) or SymPy Inventory classes. Each AbstractInventory instance
    stores nuclides and associated quantities (numbers of atoms) in the contents dictionary. A
    decay dataset is associated with each instance.

    Parameters
    ----------
    contents : dict
        Dictionary containing nuclide strings/canonical ids or Nuclide instances as keys and
        quantities (activities, number of atoms, masses or moles) as values.
    units : str, optional
        Units of the contents dictionary values. Specify either 'num' for number of atoms, or an
        activity, mass or moles unit. e.g. 'Bq', 'kBq', 'Ci', 'g', 'kg', 'mol', 'kmol'. Default is
        'Bq'.
    check : bool, optional
        Check for the validity of the contents dictionary (nuclides are in the decay dataset,
        and quantities provided are physical, etc.). Default is True.
    decay_data : DecayData, optional
        Decay dataset (default is the ICRP-107 dataset).

    Attributes
    ----------
    contents : dict
        Dictionary containing nuclide strings as keys and number of atoms of each nuclide as
        values. Nuclides are sorted alphabetically in this dictionary.
    decay_data : DecayData
        Decay dataset.
    decay_matrices : DecayMatrices
        SciPy or SymPy version of a DecayMatrices instance associated with the decay dataset.

    """

    def __init__(
        self,
        contents: Dict[Union[str, int, Nuclide], Union[float, Expr]],
        units: str = "Bq",
        check: bool = True,
        decay_data: DecayData = DEFAULTDATA,
    ) -> None:
        self.decay_data = decay_data
        self.decay_matrices = self._get_decay_matrices()

        if check is True:
            contents_with_parsed_keys: Dict[str, Union[float, Expr]] = (
                self._parse_nuclides(
                    contents, self.decay_data.nuclides, self.decay_data.dataset_name
                )
            )
            self._check_values(contents_with_parsed_keys)
        else:
            contents_with_parsed_keys = contents
        contents_sorted = sort_dictionary_alphabetically(contents_with_parsed_keys)
        self.contents = self._convert_to_number(
            contents_sorted, units, self.decay_data.dataset_name
        )

    @staticmethod
    def _parse_nuclides(
        contents: Dict[Union[str, int, Nuclide], Union[float, Expr]],
        nuclides: np.ndarray,
        dataset_name: str,
    ) -> Dict[str, Union[float, Expr]]:
        """
        Checks that nuclide keys in the contents dictionary. Converts Nuclide instances into
        nuclide name strings. Converts nuclide name strings and ids into Ab-XY format.
        """

        return {
            (
                parse_nuclide(nuc.nuclide, nuclides, dataset_name)
                if isinstance(nuc, Nuclide)
                else parse_nuclide(nuc, nuclides, dataset_name)
            ): inp
            for nuc, inp in contents.items()
        }

    @staticmethod
    def _check_values(contents: Dict[str, Union[float, Expr]]) -> None:
        """Checks that the nuclide quantities in contents (dictionary values) are valid."""

        for nuc, inp in contents.items():
            if isinstance(inp, (numbers.Number, Expr)):
                if inp >= 0:
                    continue
            raise ValueError(f"{inp} is not a valid quantity of nuclide {nuc}.")

    def _get_atomic_mass(self, nuc: str) -> Union[float, Expr]:
        """Returns the appropriate atomic mass."""

        return self.decay_matrices.atomic_masses[self.decay_data.nuclide_dict[nuc]]

    def _get_decay_const(self, nuc: str) -> Union[float, Expr]:
        """Returns the appropriate atomic mass."""

        return self.decay_matrices.decay_consts[self.decay_data.nuclide_dict[nuc]]

    @abstractmethod
    def _get_decay_matrices(self) -> DecayMatrices:
        """Returns the appropriate DecayMatrices instance."""

    @staticmethod
    @abstractmethod
    def _get_quantity_converter() -> Type[QuantityConverter]:
        """Returns the appropriate QuantityConverter instance."""

    @abstractmethod
    def _get_unit_converter(self) -> Type[UnitConverter]:
        """Returns the appropriate UnitConverter instance."""

    @abstractmethod
    def _get_year_conv(self) -> Union[float, Expr]:
        """Returns the appropriate number of days in a year."""

    def _convert_to_number(
        self,
        contents: Dict[str, Union[float, Expr]],
        units: str,
        dataset_name: str,
    ) -> Dict[str, Union[float, Expr]]:
        """
        Converts an inventory dictionary where the values are masses, moles or activities to one
        where the values are number of atoms.

        Parameters
        ----------
        contents : dict
            Dictionary containing nuclide strings or Nuclide objects as keys, and masses,
            moles, or activities as values.
        units : str
            Units of the values in the input dictionary.

        Returns
        -------
        dict
            Inventory dictionary where the values are the number of atoms of the nuclide.

        Raises
        ------
        ValueError
            If the units supplied are invalid or an activity is supplied for a stable nuclide.

        """

        if units == "num":
            contents_as_numbers = contents
        elif units in self._get_unit_converter().activity_units:
            contents_as_numbers = {}
            for nuc, act in contents.items():
                decay_const = self._get_decay_const(nuc)
                if decay_const == 0:
                    raise ValueError(
                        f"{nuc} is a stable radionuclide in {dataset_name} decay dataset. "
                        "Initializing an inventory with an activity for a stable nuclide (i.e. "
                        f"{act} {units} for {nuc}) has no physical meaning. Initialize using a "
                        "mass or mole quantity instead."
                    )
                contents_as_numbers[
                    nuc
                ] = self._get_quantity_converter().activity_to_number(
                    self._get_unit_converter().activity_unit_conv(act, units, "Bq"),
                    decay_const,
                )
        elif units in self._get_unit_converter().moles_units:
            contents_as_numbers = {
                nuc: self._get_quantity_converter().moles_to_number(
                    self._get_unit_converter().moles_unit_conv(mol, units, "mol")
                )
                for nuc, mol in contents.items()
            }
        elif units in self._get_unit_converter().mass_units:
            contents_as_numbers = {
                nuc: self._get_quantity_converter().mass_to_number(
                    self._get_unit_converter().mass_unit_conv(mass, units, "g"),
                    self._get_atomic_mass(nuc),
                )
                for nuc, mass in contents.items()
            }
        else:
            raise ValueError(f"{units} is not a supported unit.")

        return contents_as_numbers

    @property
    def nuclides(self) -> List[str]:
        """
        Returns a list of the nuclides in the inventory.

        Examples
        --------
        >>> rd.Inventory({'Tc-99m': 2.3, 'I-123': 5.8}, 'Bq').nuclides
        ['I-123', 'Tc-99m']

        """

        return list(self.contents.keys())

    def numbers(self) -> Dict[str, float]:
        """
        Returns a dictionary containing the number of atoms of each nuclide within the inventory.

        Examples
        --------
        >>> rd.Inventory({'Tc-99m': 2.3, 'I-123': 5.8}, 'Bq').numbers()
        {'I-123': 399738.47946141585, 'Tc-99m': 71852.27235544211}

        """

        return self.contents

    def activities(self, units: str = "Bq") -> Dict[str, float]:
        """
        Returns a dictionary containing the activity of each nuclide within the inventory.

        Parameters
        ----------
        units : str, optional
            Activity units for output, e.g. 'Bq', 'kBq', 'mBq', 'Ci', 'dpm'...
            Deafult is 'Bq'.

        Examples
        --------
        >>> rd.Inventory({'Tc-99m': 2300, 'I-123': 5800}, 'Bq').activities('kBq')
        {'I-123': 5.8, 'Tc-99m': 2.3}

        """

        activities = {
            nuc: self._get_unit_converter().activity_unit_conv(
                self._get_quantity_converter().number_to_activity(
                    num, self._get_decay_const(nuc)
                ),
                "Bq",
                units,
            )
            for nuc, num in self.contents.items()
        }

        return activities

    def masses(self, units: str = "g") -> Dict[str, float]:
        """
        Returns a dictionary containing the mass of each nuclide within the inventory

        Parameters
        ----------
        units : str, optional
            Mass units for output, e.g. 'Bq', 'g', 'kg', 'mg', 'ton'...
            Deafult is 'g'.

        Examples
        --------
        >>> rd.Inventory({'Tc-99m': 2.3, 'I-123': 5.8}, 'Bq').masses('pg')
        {'I-123': 8.158243973887584e-05, 'Tc-99m': 1.1800869622748502e-05}

        """

        masses = {
            nuc: self._get_unit_converter().mass_unit_conv(
                self._get_quantity_converter().number_to_mass(
                    num, self._get_atomic_mass(nuc)
                ),
                "g",
                units,
            )
            for nuc, num in self.contents.items()
        }

        return masses

    def moles(self, units: str = "mol") -> Dict[str, float]:
        """
        Returns a dictionary containing the number of atoms of each nuclide within the inventory in
        moles.

        Parameters
        ----------
        units : str, optional
            Moles units, e.g. 'mmol', 'mol', 'kmol'...
            Deafult is 'mol'.

        Examples
        --------
        >>> rd.Inventory({'Tc-99m': 2.3, 'I-123': 5.8}, 'Bq').moles()
        {'I-123': 6.637813617983513e-19, 'Tc-99m': 1.1931350531142702e-19}

        """

        moles = {
            nuc: self._get_unit_converter().moles_unit_conv(
                self._get_quantity_converter().number_to_moles(num), "mol", units
            )
            for nuc, num in self.contents.items()
        }

        return moles

    def activity_fractions(self) -> Dict[str, float]:
        """
        Returns a dictionary containing the activity fraction of each nuclide within the inventory.

        Examples
        --------
        >>> rd.Inventory({'Tc-99m': 2.3, 'I-123': 5.8}, 'Bq').activity_fractions()
        {'I-123': 0.7160493827160493, 'Tc-99m': 0.2839506172839506}

        """

        total_activity = sum(self.activities().values())
        activity_fractions = {
            nuc: act / total_activity for nuc, act in self.activities().items()
        }

        return activity_fractions

    def mass_fractions(self) -> Dict[str, float]:
        """
        Returns a dictionary containing the mass fraction of each nuclide within the inventory.

        Examples
        --------
        >>> rd.Inventory({'Tc-99m': 2.3, 'I-123': 5.8}, 'Bq').mass_fractions()
        {'I-123': 0.8736297770616593, 'Tc-99m': 0.12637022293834066}

        """

        total_mass = sum(self.masses().values())
        mass_fractions = {nuc: mass / total_mass for nuc, mass in self.masses().items()}

        return mass_fractions

    def mole_fractions(self) -> Dict[str, float]:
        """
        Returns a dictionary containing the mole fraction of each nuclide within the inventory.

        Examples
        --------
        >>> rd.Inventory({'Tc-99m': 2.3, 'I-123': 5.8}, 'Bq').mole_fractions()
        {'I-123': 0.8476385041932588, 'Tc-99m': 0.15236149580674116}

        """

        total_number = sum(self.numbers().values())
        mole_fractions = {
            nuc: num / total_number for nuc, num in self.numbers().items()
        }

        return mole_fractions

    def __len__(self) -> int:
        """
        Returns number of nuclides in the inventory.
        """

        return len(self.contents)

    def add(
        self,
        add_contents: Dict[Union[str, int, Nuclide], float],
        units: str = "Bq",
    ) -> None:
        """
        Adds a dictionary of nuclides and associated quantities (numbers/activities/masses/moles)
        to the inventory.

        Parameters
        ----------
        add_contents : dict
            Dictionary containing nuclide strings, canonical ids or Nuclide objects as keys and the
            amount of each nuclide (with specified units) as values which are added to the
            Inventory.
        units : str, optional
            Units of the values in the dictionary (e.g. 'Bq', 'Ci', 'g', 'mol', 'num'; 'Bq' is the
            default).

        Examples
        --------
        >>> inv = rd.Inventory({'H-3': 1.0}, 'Bq')
        >>> inv.add({'C-14': 2.0}, 'kBq')
        >>> inv.add({190400000: 20.0})
        >>> inv.activities()
        {'K-40': 20.0, 'C-14': 2000.0, 'H-3': 1.0}

        """

        other = Inventory(add_contents, units, True, self.decay_data)
        self.contents = (self + other).contents

    def subtract(
        self,
        sub_contents: Dict[Union[str, int, Nuclide], float],
        units: str = "Bq",
    ) -> None:
        """
        Subtracts a dictionary of nuclides and associated numbers/activities/masses
        from the inventory.

        Parameters
        ----------
        sub_contents : dict
            Dictionary containing nuclide strings or Nuclide objects as keys
            and the amount of each nuclide (with specified units) as values
            which are subtracted from the Inventory.
        units : str, optional
            Units of the values in the dictionary (e.g. 'Bq', 'Ci', 'g', 'mol', 'num'; 'Bq' is the
            default).

        Examples
        --------
        >>> inv = rd.Inventory({'K-40': 20.0, 'C-14': 2.0, 'H-3': 1.0}, 'Bq')
        >>> inv.subtract({'H-3': 1.0})
        >>> inv.subtract({190400000: 20.0})
        >>> inv.activities()
        {'C-14': 2.0, 'H-3': 0.0}

        """

        other = Inventory(sub_contents, units, True, self.decay_data)
        self.contents = (self - other).contents

    def __add__(self, other: "AbstractInventory") -> "AbstractInventory":
        """Defines + operator to add two Inventory objects together."""

        if self.decay_data != other.decay_data:
            raise ValueError(
                f"Decay datasets do not match. {self.decay_data.dataset_name} and "
                f"{other.decay_data.dataset_name}"
            )
        new_contents = add_dictionaries(self.contents, other.contents)
        return self.__class__(new_contents, "num", False, self.decay_data)

    def __sub__(self, other: "AbstractInventory") -> "AbstractInventory":
        """Defines - operator to subtract one inventory object from another."""

        if self.decay_data != other.decay_data:
            raise ValueError(
                f"Decay datasets do not match. {self.decay_data.dataset_name} and "
                f"{other.decay_data.dataset_name}"
            )
        sub_contents = other.contents.copy()
        sub_contents.update(
            (nuclide, number * -1.0) for nuclide, number in sub_contents.items()
        )
        new_contents = add_dictionaries(self.contents, sub_contents)
        return self.__class__(new_contents, "num", False, self.decay_data)

    def __mul__(self, const: Union[float, Expr]) -> "AbstractInventory":
        """
        Defines * operator to multiply all quantities of nuclides in an inventory by a constant.
        """

        new_contents = self.contents.copy()
        for nuclide, number in new_contents.items():
            new_contents[nuclide] = number * const

        return self.__class__(new_contents, "num", False, self.decay_data)

    def __rmul__(self, const: Union[float, Expr]) -> "AbstractInventory":
        """
        Defines * operator to multiply all quantities of nuclides in an Inventory by a constant.
        """

        return self.__mul__(const)

    def __truediv__(self, const: Union[float, Expr]) -> "AbstractInventory":
        """
        Defines / operator to divide all quantities of nuclides in an inventory by a constant.
        """
        new_contents = self.contents.copy()
        for nuclide, number in new_contents.items():
            new_contents[nuclide] = number / const

        return self.__class__(new_contents, "num", False, self.decay_data)

    @singledispatchmethod
    def remove(
        self, delete: Union[str, int, Nuclide, List[Union[str, int, Nuclide]]]
    ) -> None:
        """
        Removes nuclide(s) from the inventory.

        Parameters
        ----------
        delete : str or int or Nuclide or list
            Nuclide string, canonical id, Nuclide object or list of nuclide strings, canonical ids
            or Nuclide objects to delete from the Inventory object.

        Examples
        --------
        >>> inv = rd.Inventory({'Be-10': 2.0, 'C-14': 3.0, 'H-3': 1.0, 'K-40': 4.0}, 'Bq')
        >>> inv.remove('H-3')
        >>> inv.activities()
        {'Be-10': 2.0, 'C-14': 3.0, 'K-40': 4.0}
        >>> inv.remove(['Be-10', 'K-40'])
        >>> inv.activities()
        {'C-14': 3.0}

        """
        raise NotImplementedError(
            "Method takes a nuclide string, a canonical id, a Nuclide instance, or list of nuclide"
            " strings, canonical ids or Nuclide instances as a parameter."
        )

    @remove.register(str)
    def _(
        self, delete: str
    ) -> Callable[[Dict[str, float], str, bool, DecayData], None]:
        """Remove nuclide string from this inventory."""

        delete = parse_nuclide(
            delete, self.decay_data.nuclides, self.decay_data.dataset_name
        )
        new_contents = self.contents.copy()
        if delete not in new_contents:
            raise ValueError(f"{delete} does not exist in this inventory.")
        new_contents.pop(delete)

        self.contents = new_contents

    @remove.register(int)
    def _(
        self, delete: int
    ) -> Callable[[Dict[str, float], str, bool, DecayData], None]:
        """Remove nuclide string from this inventory."""

        delete_str = parse_nuclide(
            delete, self.decay_data.nuclides, self.decay_data.dataset_name
        )
        new_contents = self.contents.copy()
        if delete_str not in new_contents:
            raise ValueError(f"{delete_str} does not exist in this inventory.")
        new_contents.pop(delete_str)

        self.contents = new_contents

    @remove.register(Nuclide)
    def _(
        self, delete: Nuclide
    ) -> Callable[[Dict[str, float], str, bool, DecayData], None]:
        """Remove Nuclide object from this inventory."""

        delete_str = parse_nuclide(
            delete.nuclide, self.decay_data.nuclides, self.decay_data.dataset_name
        )
        new_contents = self.contents.copy()
        if delete_str not in new_contents:
            raise ValueError(f"{delete_str} does not exist in this inventory.")
        new_contents.pop(delete_str)

        self.contents = new_contents

    @remove.register(list)
    def _(
        self, delete: List[Union[str, int, Nuclide]]
    ) -> Callable[[Dict[str, float], str, bool, DecayData], None]:
        """Remove list of nuclide(s) from this inventory."""

        delete_list = [
            (
                parse_nuclide(
                    nuc.nuclide, self.decay_data.nuclides, self.decay_data.dataset_name
                )
                if isinstance(nuc, Nuclide)
                else parse_nuclide(
                    nuc, self.decay_data.nuclides, self.decay_data.dataset_name
                )
            )
            for nuc in delete
        ]
        new_contents = self.contents.copy()
        for nuc in delete_list:
            if nuc not in new_contents:
                raise ValueError(f"{nuc} does not exist in this inventory.")
            new_contents.pop(nuc)

        self.contents = new_contents

    def _convert_decay_time(
        self, decay_time: Union[float, Expr], units: str
    ) -> Union[float, Expr]:
        """
        Converts a decay time period into seconds.
        """

        return self._get_unit_converter().time_unit_conv(
            decay_time, units, "s", self._get_year_conv()
        )

    def _setup_decay_calc(
        self,
    ) -> Tuple[Union[np.ndarray, Matrix], List[int], Union[sparse.csr_matrix, Matrix]]:
        """
        Setup variables for a decay calculation.
        """

        vector_n0 = self.decay_matrices.vector_n0.copy()
        indices_set = set()
        for nuclide in self.contents:
            idx = self.decay_data.nuclide_dict[nuclide]
            vector_n0[idx] = self.contents[nuclide]
            indices_set.update(self.decay_data.scipy_data.matrix_c[:, idx].nonzero()[0])
        indices = list(indices_set)

        matrix_e = self.decay_matrices.matrix_e.copy()

        return vector_n0, indices, matrix_e

    def _perform_decay_calc(
        self,
        vector_n0: Union[np.ndarray, Matrix],
        matrix_e: Union[sparse.csr_matrix, Matrix],
    ) -> Union[np.ndarray, Matrix]:
        """
        Perform decay calculation matrix multiplication.
        """

        vector_nt = (
            self.decay_matrices.matrix_c
            @ matrix_e
            @ self.decay_matrices.matrix_c_inv
            @ vector_n0
        )
        return vector_nt

    @abstractmethod
    def decay(self, decay_time: float, units: str = "s") -> "AbstractInventory":
        """
        Returns a new inventory calculated from the radioactive decay of the current inventory for
        decay_time.
        """

    @abstractmethod
    def cumulative_decays(
        self, decay_time: float, units: str = "s"
    ) -> Dict[str, float]:
        """
        Calculates the total number of decays of each nuclide in the inventory between t=0 and
        t=decay_time. Note no results are reported for stable nuclides, as cumulative decays is
        zero.
        """

    def half_lives(self, units: str = "s") -> Dict[str, Union[float, str]]:
        """
        Returns dictionary of half-lives of the nuclides in the inventory in your chosen time
        units, or as a human-readable string with appropriate units.

        Parameters
        ----------
        units : str, optional
            Units for half-life. Options are 'ps', 'ns', 'μs', 'us', 'ms', 's', 'm', 'h', 'd', 'y',
            'ky', 'My', 'By', 'Gy', 'Ty', 'Py', and common spelling variations. Default is 's', i.e.
            seconds. Use 'readable' to get strings of the half-lives in human-readable units.

        Returns
        -------
        dict
            Dictionary with nuclide strings as keys and half-life floats or human-readable
            half-life strings as values.

        Examples
        --------
        >>> inv = rd.Inventory({'C-14': 1.0, 'K-40': 2.0}, 'Bq')
        >>> inv.half_lives('y')
        {'C-14': 5700.0, 'K-40': 1251000000.0}
        >>> inv.half_lives('readable')
        {'C-14': '5.70 ky', 'K-40': '1.251 By'}
        """

        return {nuc: self.decay_data.half_life(nuc, units) for nuc in self.contents}

    def progeny(self) -> Dict[str, List[str]]:
        """
        Returns dictionary with the direct progeny of the nuclides in the inventory.

        Returns
        -------
        dict
            Dictionary with nuclide strings as keys and lists of the direct progeny of each
            nuclide, ordered by decreasing branching fraction, as values.

        Examples
        --------
        >>> inv = rd.Inventory({'C-14': 1.0, 'K-40': 2.0}, 'Bq')
        >>> inv.progeny()
        {'C-14': ['N-14'], 'K-40': ['Ca-40', 'Ar-40'])

        """

        return {nuc: Nuclide(nuc, self.decay_data).progeny() for nuc in self.contents}

    def branching_fractions(self) -> Dict[str, List[float]]:
        """
        Returns dictionary with the branching fractions of the direct progeny of the nuclides
        in the inventory.

        Returns
        -------
        dict
            Dictionary with nuclide strings as keys and lists of the branching fractions to
            the direct progeny of each nuclide as values.

        Examples
        --------
        >>> inv = rd.Inventory({'C-14': 1.0, 'K-40': 2.0}, 'Bq')
        >>> inv.branching_fractions()
        {'C-14': [1.0], 'K-40': [0.8914, 0.1086])

        """

        return {
            nuc: Nuclide(nuc, self.decay_data).branching_fractions()
            for nuc in self.contents
        }

    def decay_modes(self) -> Dict[str, List[str]]:
        """
        Returns dictionary with the decay modes of the direct progeny of the nuclides in the
        inventory. Note: the decay mode strings returned are not lists of all the different
        radiation types emitted during the parent to progeny decay processes. They are the labels
        defined in the decay dataset to classify the decay type (e.g. '\u03b1', '\u03b2-' or 'IT').

        Returns
        -------
        dict
            Dictionary with nuclide strings as keys and lists of the decay modes of the parent
            nuclide

        Examples
        --------
        >>> inv = rd.Inventory({'C-14': 1.0, 'K-40': 2.0}, 'Bq')
        >>> inv.decay_modes()
        {'C-14': ['\u03b2-'], 'K-40': ['\u03b2-', '\u03b2+ \u0026 EC'])

        """

        return {
            nuc: Nuclide(nuc, self.decay_data).decay_modes() for nuc in self.contents
        }

    def decay_time_series_pandas(
        self,
        time_period: Union[float, np.ndarray],
        time_units: str = "s",
        time_scale: str = "linear",
        decay_units: str = "Bq",
        npoints: int = 501,
    ) -> pd.DataFrame:
        """
        Returns a dataframe with the initial isotope and all decay progeny decayed for the amount of time specified by
        time_period.

        Parameters
        ----------
        time_period : Union[float, np.ndarray]
            Time to decay the chain for. If a float is given, <time_scale> and <npoints> is used to create an evenly spaced array of numbers.
            If an numpy ndarray is provided, the values contained within are used and <npoints> is ignored.
        time_units : str, optional
            Units for time series. Options are 'ps', 'ns', 'μs', 'us', 'ms', 's', 'm', 'h', 'd', 'y',
            'ky', 'My', 'By', 'Gy', 'Ty', 'Py', and common spelling variations. Default is 's', i.e.
            seconds.
        time_scale : str, optional
            The time axis scale type to apply ('linear' or 'log', default is 'linear').
        decay_units : str, optional
            Units to use for the decayed values e.g. 'Bq', 'kBq', 'Ci', 'g', 'mol', 'num',
            'activity_frac', 'mass_frac', 'mol_frac'. Default is 'Bq'.
        npoints : None or int, optional
            Number of time points used to plot graph (default is 501 for normal precision decay
            calculations, or 51 for high precision decay calculations).

        Returns
        -------
        pandas.DataFrame
            Pandas DataFrame with the data of the decayed Inventory. Each isotope is its own column, with a row for
            each time increment. The time column is set as the index.

        Raises
        ------
        ValueError
            If the decay_units parameter is invalid.

        Examples
        --------
        >>> inv = rd.Inventory({'C-14': 1.0})
        >>> inv.decay_data_as_dataframe(time_period=10, time_units='ky', decay_units='mass_frac', npoints=4)
                       C-14      N-14
        Time (ky)
        0.000000   1.000000  0.000000
        3.333333   0.666747  0.333253
        6.666667   0.444550  0.555450
        10.000000  0.296402  0.703598
        """
        tmin = 0.0 if time_scale == "linear" else 0.1

        time_points = (
            time_period
            if isinstance(time_period, np.ndarray)
            else (
                np.linspace(tmin, time_period, num=npoints)
                if time_scale == "linear"
                else np.logspace(np.log10(tmin), np.log10(time_period), num=npoints)
            )
        )

        data_map = map(self.decay, time_points, itertools.repeat(time_units))
        decayed_data = collections.defaultdict(list)

        if decay_units in self._get_unit_converter().activity_units:
            for e in data_map:
                for iso, en in e.activities(decay_units).items():
                    decayed_data[iso].append(en)
        elif decay_units in self._get_unit_converter().moles_units:
            for e in data_map:
                for iso, en in e.moles(decay_units).items():
                    decayed_data[iso].append(en)
        elif decay_units in self._get_unit_converter().mass_units:
            for e in data_map:
                for iso, en in e.masses(decay_units).items():
                    decayed_data[iso].append(en)
        elif decay_units == "num":
            for e in data_map:
                for iso, en in e.numbers().items():
                    decayed_data[iso].append(en)
        elif decay_units == "activity_frac":
            for e in data_map:
                for iso, en in e.activity_fractions().items():
                    decayed_data[iso].append(en)
        elif decay_units == "mass_frac":
            for e in data_map:
                for iso, en in e.mass_fractions().items():
                    decayed_data[iso].append(en)
        elif decay_units == "mol_frac":
            for e in data_map:
                for iso, en in e.mole_fractions().items():
                    decayed_data[iso].append(en)
        else:
            raise ValueError(f"{decay_units} is not a supported decay unit.")

        time_column: str = f"Time ({time_units})"
        decayed_data[time_column] = list(time_points)

        return pd.DataFrame(decayed_data).set_index(time_column)

    def decay_time_series(
        self,
        time_period: Union[float, np.ndarray],
        time_units: str = "s",
        time_scale: str = "linear",
        decay_units: str = "Bq",
        npoints: int = 501,
    ) -> Tuple[List[float], Dict[str, List[float]]]:
        """
        Returns a list, dict tuple with the initial isotope and all decay progeny decayed for the amount of time specified by
        time_period. The list contains the time data from the decay calculations, and the dict has the isotope as the key,
        with decay data contained in a list for the value.

        Parameters
        ----------
        time_period : Union[float, np.ndarray]
            Time to decay the chain for. If a float is given, <time_scale> and <npoints> is used to create an evenly spaced array of numbers.
            If an numpy ndarray is provided, the values contained within are used and <npoints> is ignored.
        time_units : str, optional
            Units for time series. Options are 'ps', 'ns', 'μs', 'us', 'ms', 's', 'm', 'h', 'd', 'y',
            'ky', 'My', 'By', 'Gy', 'Ty', 'Py', and common spelling variations. Default is 's', i.e.
            seconds.
        time_scale : str, optional
            The time axis scale type to apply ('linear' or 'log', default is 'linear').
        decay_units : str, optional
            Units to use for the decayed values e.g. 'Bq', 'kBq', 'Ci', 'g', 'mol', 'num',
            'activity_frac', 'mass_frac', 'mol_frac'. Default is 'Bq'.
        npoints : None or int, optional
            Number of time points used to plot graph (default is 501 for normal precision decay
            calculations, or 51 for high precision decay calculations).

        Returns
        -------
        List[float]
            The time points of the decayed data
        Dict[str, List[float]]
            The isotopes as the key with the decay values as a list of floats

        Raises
        ------
        ValueError
            If the decay_units parameter is invalid.

        Examples
        --------
        >>> inv = rd.Inventory({'C-14': 1.0})
        >>> time, data = inv.decay_time_series_pandas(time_period=10, time_units="ky", decay_units="mass_frac", npoints=4)
        >>> time
        [0.0, 3.3333333333333335, 6.666666666666667, 10.0]
        >>> data
        {'C-14': [1.0, 0.6667465897861368, 0.4445504227269143, 0.2964018201597633], 'N-14': [0.0, 0.3332534102138631, 0.5554495772730857, 0.7035981798402366]}

        """
        df = self.decay_time_series_pandas(
            time_period=time_period,
            time_units=time_units,
            time_scale=time_scale,
            decay_units=decay_units,
            npoints=npoints,
        )

        return (list(df.index), df.to_dict(orient="list"))

    def plot(  # type: ignore
        self,
        xmax: float,
        xunits: str = "s",
        xmin: float = 0.0,
        xscale: str = "linear",
        yscale: str = "linear",
        ymin: float = 0.0,
        ymax: Optional[float] = None,
        yunits: str = "Bq",
        display: Union[str, List[str]] = "all",
        order: str = "dataset",
        npoints: int = 501,
        fig: Optional[matplotlib.figure.Figure] = None,
        axes: Optional[matplotlib.axes.Axes] = None,
        **kwargs,
    ) -> Tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
        """
        Plots a decay graph showing the change in activity of the inventory over time. Creates
        matplotlib fig, axes objects if they are not supplied. Returns fig, axes tuple.

        Parameters
        ----------
        xmax : float
            Maximum decay time on x-axis.
        xunits : str, optional
            Units for decay times (default is 's', i.e. seconds). Options are 'ps', 'ns', 'us',
            'ms', 's', 'm', 'h', 'd', 'y', 'ky', 'My', 'Gy', 'Ty', 'Py', and some of the common
            spelling variations of these time units.
        xmin : float, optional
            Minimum decay time on x-axis (default is 0.0 for linear x-axis, 0.1 for log x-axis).
        xscale : str, optional
            The time axis scale type to apply ('linear' or 'log', default is 'linear').
        yscale : str, optional
            The y-axis scale type to apply ('linear' or 'log', default is 'linear').
        ymin : float, optional
            Minimum value for the y-axis (default is 0.0 for ``yscale='linear'``, or 0.95x the
            minimum quantity that occurs over the decay period for ``yscale='log'``).
        ymax : None or float, optional
            Maximum value for the y-axis. Default is None, which sets the limit to 1.05x the
            maximum quantity that occurs over the decay period.
        yunits : str, optional
            Units to display on the y-axis e.g. 'Bq', 'kBq', 'Ci', 'g', 'mol', 'num',
            'activity_frac', 'mass_frac', 'mol_frac'. Default is 'Bq'.
        display : str or list, optional
            Only display the nuclides within this list on the graph. Use this parameter when
            you want to choose specific nuclide decay curves shown on the graph, either by
            supplying a string (to show one nuclide) or a list of strings (to show multiple).
            Default is 'all', which displays all nuclides present upon decay of the inventory.
        order : str, optional
            Order to display the nuclide decay curves on the graph if you do not specify the
            order via the display parameter. Default order is by 'dataset', which follows the order
            of the nuclides in the decay dataset (highest to lowest nuclides in the decay
            chains). Use 'alphabetical' if you want the nuclides to be ordered alphabetically.
        npoints : None or int, optional
            Number of time points used to plot graph (default is 501 for normal precision decay
            calculations, or 51 for high precision decay calculations).
        fig : None or matplotlib.figure.Figure, optional
            matplotlib figure object to use, or None makes ``radioactivedecay`` create one (default
            is None).
        axes : None or matplotlib.axes.Axes, optional
            matplotlib axes object to use, or None makes ``radioactivedecay`` create one (default
            is None).
        **kwargs, optional
            All additional keyword arguments to supply to matplotlib plot() function.

        Returns
        -------
        fig : matplotlib.figure.Figure
            matplotlib figure object used to plot decay chain.
        axes : matplotlib.axes.Axes
            matplotlib axes object used to plot decay chain.

        Raises
        ------
        ValueError
            If the order parameter is invalid.

        """

        if xscale == "linear":
            time_points = np.linspace(xmin, xmax, num=npoints)
        else:
            if xmin == 0.0:
                xmin = 0.1
            time_points = np.logspace(np.log10(xmin), np.log10(xmax), num=npoints)

        if display == "all":
            if order == "dataset":
                display = sort_list_according_to_dataset(
                    self.decay(0).nuclides, self.decay_data.nuclide_dict
                )
            elif order == "alphabetical":
                display = self.decay(0).nuclides
            else:
                raise ValueError(
                    f"{order} is not a valid string for the order parameter."
                )
        else:
            if isinstance(display, str):
                display = [display]
            display = [
                parse_nuclide(
                    rad, self.decay_data.nuclides, self.decay_data.dataset_name
                )
                for rad in display
            ]

        ydata = np.zeros(shape=(npoints, len(display)))
        if yunits in self._get_unit_converter().activity_units:
            for idx in range(0, npoints):
                decayed_contents = self.decay(time_points[idx], xunits).activities(
                    yunits
                )
                ydata[idx] = [decayed_contents[rad] for rad in display]
                ylabel = f"Activity ({yunits})"
        elif yunits in self._get_unit_converter().moles_units:
            for idx in range(0, npoints):
                decayed_contents = self.decay(time_points[idx], xunits).moles(yunits)
                ydata[idx] = [decayed_contents[rad] for rad in display]
                ylabel = f"Number of moles ({yunits})"
        elif yunits in self._get_unit_converter().mass_units:
            for idx in range(0, npoints):
                decayed_contents = self.decay(time_points[idx], xunits).masses(yunits)
                ydata[idx] = [decayed_contents[rad] for rad in display]
                ylabel = f"Mass ({yunits})"
        elif yunits == "num":
            for idx in range(0, npoints):
                decayed_contents = self.decay(time_points[idx], xunits).numbers()
                ydata[idx] = [decayed_contents[rad] for rad in display]
                ylabel = "Number of atoms"
        elif yunits == "activity_frac":
            for idx in range(0, npoints):
                decayed_contents = self.decay(
                    time_points[idx], xunits
                ).activity_fractions()
                ydata[idx] = [decayed_contents[rad] for rad in display]
                ylabel = "Activity fraction"
        elif yunits == "mass_frac":
            for idx in range(0, npoints):
                decayed_contents = self.decay(time_points[idx], xunits).mass_fractions()
                ydata[idx] = [decayed_contents[rad] for rad in display]
                ylabel = "Mass fraction"
        elif yunits == "mol_frac":
            for idx in range(0, npoints):
                decayed_contents = self.decay(time_points[idx], xunits).mole_fractions()
                ydata[idx] = [decayed_contents[rad] for rad in display]
                ylabel = "Mole fraction"
        else:
            raise ValueError(f"{yunits} is not a supported y-axes unit.")

        if yscale == "log" and ymin == 0.0:
            ymin = 0.95 * ydata.min()
        ylimits = [ymin, ymax] if ymax else [ymin, 1.05 * ydata.max()]

        fig, axes = decay_graph(
            time_points=time_points,
            ydata=ydata.T,
            nuclides=display,
            xunits=xunits,
            ylabel=ylabel,
            xscale=xscale,
            yscale=yscale,
            ylimits=ylimits,
            display=set(display),
            fig_in=fig,
            axes_in=axes,
            **kwargs,
        )

        return fig, axes

    def to_csv(
        self,
        filename: Union[str, pathlib.Path],
        units: str = "Bq",
        delimiter: str = ",",
        write_units: bool = False,
        header: Optional[List[str]] = None,
        encoding: str = "utf-8",
    ) -> None:
        """
        Write the contents of the inventory to a CSV file.

        The first column written is the nuclide, the second is the quantity of each nuclide in the
        chosen ``units`` (activity, mass, moles etc.), and optionally (if ``write_units=True``) the
        third column is the unit.

        Parameters
        ----------
        filename : str or pathlib.Path
            The name or path of the file to write.
        units : str, optional
            The units to output each nuclide quantity in. Default is 'Bq'. Specify 'num' if you
            require the number of atoms of each nuclide.
        delimiter : str, optional
            The delimiter to separate entries in the CSV file. Default is a comma. Specify '\t' for
            a tab-separated (i.e. TSV) file.
        write_units : bool, optional
            Write a third column with the units of each nuclide quantity. Default is False.
        header : None or str, optional
            Header row to write on first line of file. Default is None (i.e. no header row). If you
            want to include a header, pass a list of strings with the name for each column, e.g.
            ``['nuclide', 'quantity', 'units']`` if writing a units column.
        encoding : str, optional
            Encoding of the file. Default is 'utf-8'.

        Raises
        ------
        ValueError
            If the specified units are invalid.

        """

        unit_converter = self._get_unit_converter()
        if units in unit_converter.activity_units:
            contents = self.activities(units=units)
        elif units in unit_converter.mass_units:
            contents = self.masses(units=units)
        elif units in unit_converter.moles_units:
            contents = self.moles(units=units)
        elif units == "num":
            contents = self.numbers()
        else:
            raise ValueError(f"Unknown units string specified: {units}")

        rows = [header] if header else []
        for nuclide, quantity in contents.items():
            rows.append(
                [nuclide, str(quantity), units]
                if write_units
                else [nuclide, str(quantity)]
            )

        _write_csv_file(filename, rows, delimiter, encoding)

    def __eq__(self, other: object) -> bool:
        """
        Check whether two instances are equal with ``==`` operator.
        """

        if not isinstance(other, AbstractInventory):
            return NotImplemented
        return self.contents == other.contents and self.decay_data == other.decay_data

    def __ne__(self, other: object) -> bool:
        """
        Check whether two instances are not equal with ``!=`` operator.
        """

        if not isinstance(other, AbstractInventory):
            return NotImplemented
        return not self.__eq__(other)


class Inventory(AbstractInventory):
    # pylint: disable=line-too-long
    """
    ``Inventory`` instances store nuclides and associated quantities (numbers of atoms) in the
    contents dictionary. A decay dataset is associated with each instance.

    Parameters
    ----------
    contents : dict
        Dictionary containing nuclide strings/canonical ids or Nuclide instances as keys and
        quantities (activities, number of atoms, masses or moles) as values.
    units : str, optional
        Units of the contents dictionary values. Specify either 'num' for number of atoms, or an
        activity, mass or moles unit. e.g. 'Bq', 'kBq', 'Ci', 'g', 'kg', 'mol', 'kmol'. Default is
        'Bq'.
    check : bool, optional
        Check for the validity of the contents dictionary (nuclides are in the decay dataset,
        and quantities provided are physical, etc.). Default is True.
    decay_data : DecayData, optional
        Decay dataset (default is the ICRP-107 dataset).

    Attributes
    ----------
    contents : dict
        Dictionary containing nuclide strings as keys and number of atoms of each nuclide as
        values. Nuclides are sorted alphabetically in this dictionary.
    decay_data : DecayData
        Decay dataset.
    decay_matrices : DecayMatricesScipy
        Float/SciPy version of the DecayMatrices associated with the decay dataset.

    Examples
    --------
    >>> rd.Inventory({'Tc-99m': 2.3, 'I-123': 5.8}, 'Bq')
    Inventory activities (Bq): {'I-123': 2.3, 'Tc-99m': 5.8}, decay dataset: icrp107_ame2020_nubase2020
    >>> H3 = rd.Nuclide('H-3')
    >>> rd.Inventory({H3: 3.0}, 'g')
    Inventory activities (Bq): {'H-3': 1067957043281807.0}, decay dataset: icrp107_ame2020_nubase2020
    >>> rd.Inventory({270570000: 7.2, 922380000: 21.1}, 'Ci')
    Inventory activities (Bq): {'Co-57': 266400000000.0, 'U-238': 780700000000.0001}, decay dataset: icrp107_ame2020_nubase2020

    """
    # pylint: enable=line-too-long

    def _get_decay_matrices(self) -> DecayMatricesScipy:
        """Returns the appropriate DecayMatrices instance."""

        return self.decay_data.scipy_data

    @staticmethod
    def _get_quantity_converter() -> Type[QuantityConverter]:
        """
        Returns the appropriate QuantityConverter instance.
        """

        return QuantityConverterFloat

    def _get_unit_converter(self) -> Type[UnitConverter]:
        """
        Returns the appropriate UnitConverter instance.
        """

        return UnitConverterFloat

    def _get_year_conv(self) -> float:
        """Returns the appropriate number of days in a year."""

        return self.decay_data.float_year_conv

    def decay(self, decay_time: float, units: str = "s") -> "Inventory":
        """
        Returns a new inventory calculated from the radioactive decay of the current inventory for
        decay_time.

        Parameters
        ----------
        decay_time : float
            Decay time.
        units : str, optional
            Units of decay_time (default is seconds). Options are 'ns', 'us', 'ms', 's', 'm', 'h',
            'd', 'y', 'ky', 'My', 'Gy', 'Ty', 'Py', and some of the common spelling variations of
            these time units.

        Returns
        -------
        Inventory
            New Inventory after the radioactive decay.

        Examples
        --------
        >>> inv_t0 = rd.Inventory({'H-3': 1.0}, 'Bq')
        >>> inv_t1 = inv_t0.decay(12.32, 'y')
        >>> inv_t1.activities()
        {'H-3': 0.5, 'He-3': 0.0}

        """

        decay_time = self._convert_decay_time(decay_time, units)
        vector_n0, indices, matrix_e = self._setup_decay_calc()

        matrix_e.data[indices] = np.exp(
            -decay_time * self.decay_matrices.decay_consts[indices]
        )

        vector_nt = self._perform_decay_calc(vector_n0, matrix_e)
        new_contents = dict(zip(self.decay_data.nuclides[indices], vector_nt[indices]))

        return self.__class__(new_contents, "num", False, self.decay_data)

    def cumulative_decays(
        self, decay_time: float, units: str = "s"
    ) -> Dict[str, float]:
        """
        Calculates the total number of decays of each nuclide in the inventory between t=0 and
        t=decay_time. Note no results are reported for stable nuclides, as cumulative decays is
        zero.

        Parameters
        ----------
        decay_time : float
            Decay time (calculates total number of decays over this period).
        units : str, optional
            Units of decay_time (default is seconds). Options are 'ns', 'us', 'ms', 's', 'm', 'h',
            'd', 'y', 'ky', 'My', 'Gy', 'Ty', 'Py', and some of the common spelling variations of
            these time units.

        Returns
        -------
        dict
            Dictionary containing nuclide strings as keys and total number of decays of each
            nuclide as values (floats).

        Examples
        --------
        >>> inv_t0 = rd.Inventory({'Sr-90': 10.0}, 'num')
        >>> inv_t0.cumulative_decays(1.0, 'My')
        {'Sr-90': 10.0, 'Y-90': 10.000000000000002}

        """

        decay_time = self._convert_decay_time(decay_time, units)
        vector_n0, indices, matrix_e = self._setup_decay_calc()

        indices = [
            idx for idx in indices if self.decay_matrices.decay_consts[idx] > 0.0
        ]
        for idx in indices:
            matrix_e[idx, idx] = (
                1.0 - np.exp((-decay_time * self.decay_matrices.decay_consts[idx]))
            ) / self.decay_matrices.decay_consts[idx]

        cumulative_decays = self._perform_decay_calc(vector_n0, matrix_e)
        result_dict = {
            self.decay_data.nuclides[idx]: float(
                self.decay_matrices.decay_consts[idx] * cumulative_decays[idx]
            )
            for idx in indices
        }

        return result_dict

    def __repr__(self) -> str:
        return (
            f"Inventory activities (Bq): {self.activities()}, "
            f"decay dataset: {self.decay_data.dataset_name}"
        )


class InventoryHP(AbstractInventory):
    """
    ``InventoryHP`` instances store a dictionary of nuclides and associated numbers of atoms, and a
    ``DecayData`` instance of radioactive decay data. Uses SymPy high precision arithmetic for all
    calculations.

    Parameters
    ----------
    contents : dict
        Dictionary containing nuclide strings/canonical ids or Nuclide
        objects as keys and quantities as values.
    units : str, optional
        Units of the values in the contents dictionary e.g. 'Bq', 'kBq', 'Ci', 'g', 'mol',
        'num'... (default is 'Bq').
    check : bool, optional
        Check for the validity of contents and that the supplied decay dataset contains SymPy data
        (default is True).
    data : DecayData, optional
        Decay dataset (default is the ICRP-107 dataset).
    sympy_contents : dict, optional
        Version of the contents dictionary with SymPy expressions as values. Setting this requires
        that data (the decay dataset) contains SymPy data. ``radioactivedecay`` will create
        sympy_contents automatically from contents if None is supplied.

    Attributes
    ----------
    contents : dict
        Dictionary containing nuclide strings as keys and number of atoms of each nuclide as
        values. Nuclides are sorted alphabetically in this dictionary.
    data : DecayData
        Decay dataset.
    decay_matrices : DecayMatricesSympy
       SymPy DecayMatrices instance associated with the decay dataset.
    sig_fig: int
        Number of significant figures for high precision decay calculations and plots. Deafult is
        320.
    unit_converter : UnitConverterSympy
        SymPy version of a convertor for within different units.

    Examples
    --------
    >>> rd.InventoryHP({'Tc-99m': 2.3, 'I-123': 5.8}, 'Bq')
    InventoryHP activities: {'I-123': 2.3, 'Tc-99m': 5.8}, decay dataset: icrp107
    >>> H3 = rd.Nuclide('H-3')
    >>> rd.InventoryHP({H3: 3.0} 'Bq')
    InventoryHP activities: {'H-3': 3.0}, decay dataset: icrp107
    >>> rd.InventoryHP({'U-238': 21.1, 'Co-57': 7.2}, 'Ci')
    InventoryHP activities: {'Co-57': 7.2, 'U-238': 21.1}, decay dataset: icrp107

    """

    def __init__(
        self,
        contents: Dict[Union[str, int, Nuclide], float],
        units: str = "Bq",
        check: bool = True,
        decay_data: DecayData = DEFAULTDATA,
    ) -> None:
        if check is True:
            try:
                assert decay_data.sympy_data
                assert decay_data.sympy_year_conv
            except ValueError:
                raise ValueError(
                    f"Decay dataset supplied to {self.__class__.__name__} constructor does not "
                    "contain SymPy data."
                ) from None
            contents = {nuc: nsimplify(val) for nuc, val in contents.items()}

        self.sig_fig = 320
        super().__init__(contents, units, check, decay_data)

    def _get_decay_matrices(self) -> DecayMatricesSympy:
        """
        Returns the appropriate DecayMatrices instance.
        """

        return self.decay_data.sympy_data

    @staticmethod
    def _get_quantity_converter() -> Type[QuantityConverter]:
        """
        Returns the appropriate QuantityConverter instance.
        """

        return QuantityConverterSympy

    def _get_unit_converter(self) -> Type[UnitConverter]:
        """
        Returns the appropriate UnitConverter instance.
        """

        return UnitConverterSympy

    def _get_year_conv(self) -> Expr:
        """Returns the appropriate number of days in a year."""

        return self.decay_data.sympy_year_conv

    def numbers(self) -> Dict[str, float]:
        """
        Returns a dictionary containing the number of atoms of each nuclide (as floats) within this
        InventoryHP instance.

        Examples
        --------
        >>> rd.InventoryHP({'Tc-99m': 2.3, 'I-123': 5.8}, 'Bq').numbers()
        {'I-123': 399738.47946141585, 'Tc-99m': 71852.27235544211}

        """

        contents_in_floats = {nuc: float(num) for nuc, num in self.contents.items()}
        return contents_in_floats

    def activities(self, units: str = "Bq") -> Dict[str, float]:
        """
        Returns a dictionary containing the activity of each nuclide (as floats) within this
        InventoryHP instance.

        Parameters
        ----------
        units : str, optional
            Activity units for output, e.g. 'Bq', 'kBq', 'mBq', 'Ci', 'dpm'...
            Deafult is 'Bq'.

        """

        activities = {
            nuc: float(
                self._get_unit_converter().activity_unit_conv(
                    self._get_quantity_converter().number_to_activity(
                        num, self._get_decay_const(nuc)
                    ),
                    "Bq",
                    units,
                )
            )
            for nuc, num in self.contents.items()
        }

        return activities

    def masses(self, units: str = "g") -> Dict[str, float]:
        """
        Returns a dictionary containing the mass of each nuclide (as floats) within this InventoryHP
        instance.

        Parameters
        ----------
        units : str, optional
            Mass units for output, e.g. 'Bq', 'g', 'kg', 'mg', 'ton'...
            Deafult is 'g'.

        """

        masses = {
            nuc: float(
                self._get_unit_converter().mass_unit_conv(
                    self._get_quantity_converter().number_to_mass(
                        num, self._get_atomic_mass(nuc)
                    ),
                    "g",
                    units,
                )
            )
            for nuc, num in self.contents.items()
        }

        return masses

    def moles(self, units: str = "mol") -> Dict[str, float]:
        """
        Returns a dictionary containing the number of atoms of each nuclide (as floats) within this
        InventoryHP instance.

        Parameters
        ----------
        units : str, optional
            Moles units, e.g. 'mmol', 'mol', 'kmol'...
            Deafult is 'mol'.

        """

        moles = {
            nuc: float(
                self._get_unit_converter().moles_unit_conv(
                    self._get_quantity_converter().number_to_moles(num), "mol", units
                )
            )
            for nuc, num in self.contents.items()
        }

        return moles

    def decay(self, decay_time: float, units: str = "s") -> "InventoryHP":
        """
        Decay calculation with high numerical precision. This uses SymPy arbitrary-precision
        arithmetic functions for the decay calculation. The results can be more accurate than
        a normal double precision float calculation with the ``Inventory`` class when the decay
        chains contain radionuclides with very similar half-lives or half-lives that differ by many
        orders of magnitude.

        Parameters
        ----------
        decay_time : float
            Decay time.
        units : str, optional
            Units of decay_time (default is seconds). Options are 'ns', 'us', 'ms', 's', 'm', 'h',
            'd', 'y', 'ky', 'My', 'Gy', 'Ty', 'Py', and some of the common spelling variations of
            these time units.

        Returns
        -------
        InventoryHP
            New high precision inventory after the radioactive decay.

        Raises
        ------
        ValueError
            If self.sig_fig is set to be lower than 1.

        Examples
        --------
        >>> inv_t0 = rd.InventoryHP({'Fm-257': 1.0})
        >>> inv_t1 = inv_t0.decay(10.0, 'd')
        >>> inv_t1.activities()
        {'Ac-225': 0.0,
        'Am-241': 9.985270042416324e-24,
        'Am-245': 5.4061684195880344e-09,
        ...
        'Fm-257': 0.9333548028364793,
        ...
        }

        """

        if self.sig_fig < 1:
            raise ValueError(
                "InventoryHP.sig_fig attribute needs to be int greater than 0."
            )

        decay_time = nsimplify(decay_time)
        decay_time = self._convert_decay_time(decay_time, units)
        vector_n0, indices, matrix_e = self._setup_decay_calc()

        for idx in indices:
            matrix_e[idx, idx] = exp(
                (-decay_time * self.decay_matrices.decay_consts[idx]).evalf(
                    self.sig_fig
                )
            )

        vector_nt = self._perform_decay_calc(vector_n0, matrix_e)
        new_contents = {
            self.decay_data.nuclides[idx]: vector_nt[idx] for idx in indices
        }

        return self.__class__(new_contents, "num", False, self.decay_data)

    def cumulative_decays(
        self, decay_time: float, units: str = "s"
    ) -> Dict[str, float]:
        """
        Calculates the total number of decays of each nuclide in the inventory between t=0 and
        t=decay_time. Uses SymPy high precision calculations. Note no results are reported for
        stable nuclides, as cumulative decays is zero.

        Parameters
        ----------
        decay_time : float
            Decay time (calculates total number of decays over this period).
        units : str, optional
            Units of decay_time (default is seconds). Options are 'ns', 'us', 'ms', 's', 'm', 'h',
            'd', 'y', 'ky', 'My', 'Gy', 'Ty', 'Py', and some of the common spelling variations of
            these time units.

        Returns
        -------
        dict
            Dictionary containing nuclide strings as keys and total number of decays of each
            nuclide as values (floats).

        Raises
        ------
        ValueError
            If self.sig_fig is set to be lower than 1.

        Examples
        --------
        >>> inv_t0 = rd.Inventory({'Sr-90': 10.0}, 'num')
        >>> inv_t0.cumulative_decays(1.0, 'My')
        {'Sr-90': 10.0, 'Y-90': 10.0}

        """

        if self.sig_fig < 1:
            raise ValueError(
                "InventoryHP.sig_fig attribute needs to be int greater than 0."
            )

        decay_time = nsimplify(decay_time)
        decay_time = self._convert_decay_time(decay_time, units)
        vector_n0, indices, matrix_e = self._setup_decay_calc()

        indices = [
            idx for idx in indices if self.decay_matrices.decay_consts[idx] > Integer(0)
        ]
        for idx in indices:
            matrix_e[idx, idx] = (
                Integer(1)
                - exp(
                    (-decay_time * self.decay_matrices.decay_consts[idx]).evalf(
                        self.sig_fig
                    )
                )
            ) / self.decay_matrices.decay_consts[idx]

        cumulative_decays = self._perform_decay_calc(vector_n0, matrix_e)
        result_dict = {
            self.decay_data.nuclides[idx]: float(
                self.decay_matrices.decay_consts[idx] * cumulative_decays[idx]
            )
            for idx in indices
        }

        return result_dict

    # Redefines one default parameter of the plot function (npoints) for performance reasons.
    # Note rewriting all the parameters like this violates the Don't Repeat Yourself principle, but
    # currently no better solution while still maintaining auto-complete functionality for
    # the method signature with some editors.
    def plot(  # type: ignore
        self,
        xmax: float,
        xunits: str = "s",
        xmin: float = 0.0,
        xscale: str = "linear",
        yscale: str = "linear",
        ymin: float = 0.0,
        ymax: Optional[float] = None,
        yunits: str = "Bq",
        display: Union[str, List[str]] = "all",
        order: str = "dataset",
        npoints: int = 51,
        fig: Optional[matplotlib.figure.Figure] = None,
        axes: Optional[matplotlib.axes.Axes] = None,
        **kwargs,
    ) -> Tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
        """
        Plots a decay graph using high precision decay calculations. Only difference from normal
        precision decay plot (``Inventory.plot()``) is default npoints=51.

        Parameters
        ----------
        xmax : float
            Maximum decay time on x-axis.
        xunits : str, optional
            Units for decay times (default is 's', i.e. seconds). Options are 'ps', 'ns', 'us',
            'ms', 's', 'm', 'h', 'd', 'y', 'ky', 'My', 'Gy', 'Ty', 'Py', and some of the common
            spelling variations of these time units.
        xmin : float, optional
            Minimum decay time on x-axis (default is 0.0 for linear x-axis, 0.1 for log x-axis).
        xscale : str, optional
            The time axis scale type to apply ('linear' or 'log', default is 'linear').
        yscale : str, optional
            The y-axis scale type to apply ('linear' or 'log', default is 'linear').
        ymin : float, optional
            Minimum value for the y-axis (default is 0.0 for ``yscale='linear'``, or 0.95x the
            minimum quantity that occurs over the decay period for ``yscale='log'``).
        ymax : None or float, optional
            Maximum value for the y-axis. Default is None, which sets the limit to 1.05x the
            maximum quantity that occurs over the decay period.
        yunits : str, optional
            Units to display on the y-axis e.g. 'Bq', 'kBq', 'Ci', 'g', 'mol', 'num',
            'activity_frac', 'mass_frac', 'mol_frac'. Default is 'Bq'.
        display : str or list, optional
            Only display the nuclides within this list on the graph. Use this parameter when
            you want to choose specific nuclide decay curves shown on the graph, either by
            supplying a string (to show one nuclide) or a list of strings (to show multiple).
            Default is 'all', which displays all nuclides present upon decay of the inventory.
        order : str, optional
            Order to display the nuclide decay curves on the graph if you do not specify the
            order via the display parameter. Default order is by 'dataset', which follows the order
            of the nuclides in the decay dataset (highest to lowest nuclides in the decay
            chains). Use 'alphabetical' if you want the nuclides to be ordered alphabetically.
        npoints : None or int, optional
            Number of time points used to plot graph. Default is 51.
        fig : None or matplotlib.figure.Figure, optional
            matplotlib figure object to use, or None makes ``radioactivedecay`` create one (default
            is None).
        axes : None or matplotlib.axes.Axes, optional
            matplotlib axes object to use, or None makes ``radioactivedecay`` create one (default
            is None).
        **kwargs, optional
            All additional keyword arguments to supply to matplotlib plot() function.

        Returns
        -------
        fig : matplotlib.figure.Figure
            matplotlib figure object used to plot decay chain.
        axes : matplotlib.axes.Axes
            matplotlib axes object used to plot decay chain.

        Raises
        ------
        ValueError
            If the order parameter is invalid.

        """

        return super().plot(
            xmax,
            xunits,
            xmin,
            xscale,
            yscale,
            ymin,
            ymax,
            yunits,
            display,
            order,
            npoints,
            fig,
            axes,
            **kwargs,
        )

    def __repr__(self) -> str:
        return (
            f"InventoryHP activities (Bq): {self.activities()}, "
            f"decay dataset: {self.decay_data.dataset_name}"
        )
