"""
The converters module contains classes for unit and quantity conversions.

The docstring code examples assume that ``radioactivedecay`` has been imported
as ``rd``:

.. highlight:: python
.. code-block:: python

    >>> import radioactivedecay as rd

"""

from abc import ABC
from typing import Dict, Union

from sympy import Integer, nsimplify
from sympy.core.expr import Expr

AVOGADRO = 6.02214076e23


class UnitConverter(ABC):
    """
    Template class for unit converters using either floats or SymPy arithmetic.


    Class Attributes
    ----------
    time_units : dict
        Dictionary containing numbers of each time unit per second.
    activity_units : dict
        Dictionary containing amounts of each activity unit per Bq.
    mass_units : dict
        Dictionary containing amounts of each mass unit per g.
    moles_units : dict
        Dictionary containing amounts of each mole unit per mol.
    year_units : set
        Set containing all year units. Used to pick up where days in year conversion is needed.
    time_unit_err_msg : str
        Error message for time unit ValueError.
    activity_unit_err_msg : str
        Error message for activity unit ValueError.
    mass_unit_err_msg : str
        Error message for mass unit ValueError.
    moles_unit_err_msg : str
        Error message for moles unit ValueError.

    """

    time_units: Dict[str, Union[float, Expr]]
    activity_units: Dict[str, Union[float, Expr]]
    mass_units: Dict[str, Union[float, Expr]]
    moles_units: Dict[str, Union[float, Expr]]
    year_units = {"y", "yr", "year", "years", "ky", "My", "By", "Gy", "Ty", "Py"}

    time_unit_err_msg = 'is not a valid time unit, e.g. "s", "m", "h", "d" or "y".'
    activity_unit_err_msg = 'is not a valid activitiy unit, e.g. "Bq", "kBq", "Ci"...'
    mass_unit_err_msg = 'is not a valid mass unit, e.g. "g", "kg", "mg"...'
    moles_unit_err_msg = 'is not a valid moles unit, e.g. "mol", "kmol", "mmol"...'

    @classmethod
    def time_unit_conv(
        cls,
        time_period: Union[float, Expr],
        units_from: str,
        units_to: str,
        year_conv: Union[float, Expr],
    ) -> Union[float, Expr]:
        """
        Converts a time period from one time unit to another.

        Parameters
        ----------
        time_period : float or Expr
            Time period before conversion.
        units_from : str
            Time unit before conversion.
        units_to : str
            Time unit after conversion.
        year_conv : Union[float, Expr]
            Number of days in one year.

        Returns
        -------
        float or Expr
            Time period in new units.

        Raises
        ------
        ValueError
            If one of the time unit parameters is invalid.

        """

        if units_from not in cls.time_units:
            raise ValueError(f"{units_from} {cls.time_unit_err_msg}")
        if units_to not in cls.time_units:
            raise ValueError(f"{units_to} {cls.time_unit_err_msg}")

        factor_from = cls.time_units[units_from]
        factor_to = cls.time_units[units_to]
        if units_from in cls.year_units:
            factor_from *= year_conv
        if units_to in cls.year_units:
            factor_to *= year_conv
        return time_period * factor_from / factor_to

    @classmethod
    def activity_unit_conv(
        cls, activity: Union[float, Expr], units_from: str, units_to: str
    ) -> Union[float, Expr]:
        """
        Converts an activity from one unit to another.

        Parameters
        ----------
        activity : float or Expr
            Activity before conversion.
        units_from : str
            Activity unit before conversion.
        units_to : str
            Activity unit after conversion.

        Returns
        -------
        float or Expr
            Activity in new units.

        Raises
        ------
        ValueError
            If one of the activity unit parameters is invalid.

        """

        if units_from not in cls.activity_units:
            raise ValueError(f"{units_from} {cls.activity_unit_err_msg}")
        if units_to not in cls.activity_units:
            raise ValueError(f"{units_to} {cls.activity_unit_err_msg}")

        return activity * cls.activity_units[units_from] / cls.activity_units[units_to]

    @classmethod
    def mass_unit_conv(
        cls, mass: Union[float, Expr], units_from: str, units_to: str
    ) -> Union[float, Expr]:
        """
        Converts a mass from one unit to another.

        Parameters
        ----------
        mass : float or Expr
            Mass before conversion.
        units_from : str
            Mass unit before conversion.
        units_to : str
            Mass unit after conversion.

        Returns
        -------
        float or Expr
            Mass in new units.

        Raises
        ------
        ValueError
            If one of the mass unit parameters is invalid.

        """

        if units_from not in cls.mass_units:
            raise ValueError(f"{units_from} {cls.mass_unit_err_msg}")
        if units_to not in cls.mass_units:
            raise ValueError(f"{units_to} {cls.mass_unit_err_msg}")

        return mass * cls.mass_units[units_from] / cls.mass_units[units_to]

    @classmethod
    def moles_unit_conv(
        cls, moles: Union[float, Expr], units_from: str, units_to: str
    ) -> Union[float, Expr]:
        """
        Converts a number of moles from one order of magnitude to another.

        Parameters
        ----------
        moles : float or Expr
            Amount before conversion.
        units_from : str
            Unit before conversion.
        units_to : str
            Unit after conversion.

        Returns
        -------
        float or Expr
            Result in chosen units.

        Raises
        ------
        ValueError
            If one of the unit parameters is invalid.

        """

        if units_from not in cls.moles_units:
            raise ValueError(f"{units_from} {cls.moles_unit_err_msg}")
        if units_to not in cls.moles_units:
            raise ValueError(f"{units_to} {cls.moles_unit_err_msg}")

        return moles * cls.moles_units[units_from] / cls.moles_units[units_to]


class UnitConverterFloat(UnitConverter):
    """
    Unit converter using floats.
    """

    time_units: Dict[str, float] = {
        "ps": 1.0e-12,
        "ns": 1.0e-9,
        "μs": 1.0e-6,
        "us": 1.0e-6,
        "ms": 1.0e-3,
        "s": 1.0,
        "m": 60.0,
        "h": 3600.0,
        "d": 86400.0,
        "y": 86400.0,
        "sec": 1,
        "second": 1,
        "seconds": 1,
        "hr": 3600.0,
        "hour": 3600.0,
        "hours": 3600.0,
        "day": 86400.0,
        "days": 86400.0,
        "yr": 86400.0,
        "year": 86400.0,
        "years": 86400.0,
        "ky": 86400.0 * 1.0e3,
        "My": 86400.0 * 1.0e6,
        "By": 86400.0 * 1.0e9,
        "Gy": 86400.0 * 1.0e9,
        "Ty": 86400.0 * 1.0e12,
        "Py": 86400.0 * 1.0e15,
    }

    activity_units: Dict[str, float] = {
        "pBq": 1.0e-12,
        "nBq": 1.0e-9,
        "μBq": 1.0e-6,
        "uBq": 1.0e-6,
        "mBq": 1.0e-3,
        "Bq": 1.0,
        "kBq": 1.0e3,
        "MBq": 1.0e6,
        "GBq": 1.0e9,
        "TBq": 1.0e12,
        "PBq": 1.0e15,
        "EBq": 1.0e18,
        "pCi": 1.0e-12 * 3.7e10,
        "nCi": 1.0e-9 * 3.7e10,
        "μCi": 1.0e-6 * 3.7e10,
        "uCi": 1.0e-6 * 3.7e10,
        "mCi": 1.0e-3 * 3.7e10,
        "Ci": 1.0 * 3.7e10,
        "kCi": 1.0e3 * 3.7e10,
        "MCi": 1.0e6 * 3.7e10,
        "GCi": 1.0e9 * 3.7e10,
        "TCi": 1.0e12 * 3.7e10,
        "PCi": 1.0e15 * 3.7e10,
        "ECi": 1.0e18 * 3.7e10,
        "dpm": 1.0 / 60.0,
    }

    mass_units: Dict[str, float] = {
        "pg": 1.0e-12,
        "ng": 1.0e-9,
        "μg": 1.0e-6,
        "ug": 1.0e-6,
        "mg": 1.0e-3,
        "g": 1.0,
        "kg": 1.0e3,
        "Mg": 1.0e6,
        "t": 1.0e6,
        "ton": 1.0e6,
    }

    moles_units: Dict[str, float] = {
        "pmol": 1.0e-12,
        "nmol": 1.0e-9,
        "μmol": 1.0e-6,
        "umol": 1.0e-6,
        "mmol": 1.0e-3,
        "mol": 1.0,
        "kmol": 1.0e3,
        "Mmol": 1.0e6,
    }

    @classmethod
    def __repr__(cls) -> str:
        return f"{cls.__name__} using double-precision floats."


class UnitConverterSympy(UnitConverter):
    """
    Unit converter using SymPy arbitrary precision operations.
    """

    time_units: Dict[str, Expr] = {
        "ps": Integer(1) / 1000000000000,
        "ns": Integer(1) / 1000000000,
        "μs": Integer(1) / 1000000,
        "us": Integer(1) / 1000000,
        "ms": Integer(1) / 1000,
        "s": Integer(1),
        "m": Integer(60),
        "h": Integer(3600),
        "d": Integer(86400),
        "y": Integer(86400),
        "sec": Integer(1),
        "second": Integer(1),
        "seconds": Integer(1),
        "hr": Integer(3600),
        "hour": Integer(3600),
        "hours": Integer(3600),
        "day": Integer(86400),
        "days": Integer(86400),
        "yr": Integer(86400),
        "year": Integer(86400),
        "years": Integer(86400),
        "ky": Integer(86400) * 1000,
        "My": Integer(86400) * 1000000,
        "By": Integer(86400) * 1000000000,
        "Gy": Integer(86400) * 1000000000,
        "Ty": Integer(86400) * 1000000000000,
        "Py": Integer(86400) * 1000000000000000,
    }

    activity_units: Dict[str, Expr] = {
        "pBq": Integer(1) / 1000000000000,
        "nBq": Integer(1) / 1000000000,
        "μBq": Integer(1) / 1000000,
        "uBq": Integer(1) / 1000000,
        "mBq": Integer(1) / 1000,
        "Bq": Integer(1),
        "kBq": Integer(1000),
        "MBq": Integer(1000000),
        "GBq": Integer(1000000000),
        "TBq": Integer(1000000000000),
        "PBq": Integer(1000000000000000),
        "EBq": Integer(1000000000000000000),
        "pCi": Integer(1) / 1000000000000 * 37000000000,
        "nCi": Integer(1) / 1000000000 * 37000000000,
        "μCi": Integer(1) / 1000000 * 37000000000,
        "uCi": Integer(1) / 1000000 * 37000000000,
        "mCi": Integer(1) / 1000 * 37000000000,
        "Ci": Integer(1) * 37000000000,
        "kCi": Integer(1000) * 37000000000,
        "MCi": Integer(1000000) * 37000000000,
        "GCi": Integer(1000000000) * 37000000000,
        "TCi": Integer(1000000000000) * 37000000000,
        "PCi": Integer(1000000000000000) * 37000000000,
        "ECi": Integer(1000000000000000000) * 37000000000,
        "dpm": Integer(1) / 60,
    }

    mass_units: Dict[str, Expr] = {
        "pg": Integer(1) / 1000000000000,
        "ng": Integer(1) / 1000000000,
        "μg": Integer(1) / 1000000,
        "ug": Integer(1) / 1000000,
        "mg": Integer(1) / 1000,
        "g": Integer(1),
        "kg": Integer(1000),
        "Mg": Integer(1000000),
        "t": Integer(1000000),
        "ton": Integer(1000000),
    }

    moles_units: Dict[str, Expr] = {
        "pmol": Integer(1) / 1000000000000,
        "nmol": Integer(1) / 1000000000,
        "μmol": Integer(1) / 1000000,
        "umol": Integer(1) / 1000000,
        "mmol": Integer(1) / 1000,
        "mol": Integer(1),
        "kmol": Integer(1000),
        "Mmol": Integer(1000000),
    }

    @classmethod
    def __repr__(cls) -> str:
        return f"{cls.__name__} using SymPy arbitrary precision calculations."


class QuantityConverter(ABC):
    """
    Template class for quantity converters using either floats or SymPy arithmetic. Converts
    activity in Bq, mass in g and moles in mol to number of atoms, and vice versa.

    Class Attributes
    ----------
    avogadro : float or Expr
        Avogadro constant (number of atoms/mol).

    """

    avogadro: Union[float, Expr]

    @staticmethod
    def activity_to_number(
        activity: Union[float, Expr], decay_const: Union[float, Expr]
    ) -> Union[float, Expr]:
        """
        Converts an activity in Bq to the number of atoms.

        Parameters
        ----------
        activity : float or sympy.core.expr.Expr
            The activity in Bq of the nuclide to be converted.
        decay_const : float or sympy.core.expr.Expr
            Decay constant of the nuclide in s^-1.

        Returns
        -------
        float or sympy.core.expr.Expr
            Number of atoms of the nuclide.

        """

        return activity / decay_const

    @classmethod
    def mass_to_number(
        cls, mass: Union[float, Expr], atomic_mass: Union[float, Expr]
    ) -> Union[float, Expr]:
        """
        Converts a mass in grams to number of atoms.

        Parameters
        ----------
        mass : float or sympy.core.expr.Expr
            The mass of the nuclide to be converted in grams.
        atomic_mass : float or sympy.core.expr.Expr
            Atomic mass of the nuclide in g/mol.

        Returns
        -------
        float or sympy.core.expr.Expr
            Number of atoms of the nuclide.

        """

        return mass / atomic_mass * cls.avogadro

    @classmethod
    def moles_to_number(cls, moles: Union[float, Expr]) -> Union[float, Expr]:
        """
        Converts number of moles to number of atoms.

        Parameters
        ----------
        moles : float or sympy.core.expr.Expr
            Number of moles to be converted.

        Returns
        -------
        float or sympy.core.expr.Expr
            Number of atoms.

        """

        return moles * cls.avogadro

    @staticmethod
    def number_to_activity(
        number: Union[float, Expr], decay_const: Union[float, Expr]
    ) -> Union[float, Expr]:
        """
        Converts number of atoms to activity in Bq.

        Parameters
        ----------
        number : float or sympy.core.expr.Expr
            The number of atoms of nuclide to be converted.
        decay_const : float or sympy.core.expr.Expr
            Decay constant of the nuclide in s^-1.

        Returns
        -------
        float or sympy.core.expr.Expr
            Activity of the nuclide in Bq.

        """

        return number * decay_const

    @classmethod
    def number_to_mass(
        cls, number: Union[float, Expr], atomic_mass: Union[float, Expr]
    ) -> Union[float, Expr]:
        """
        Converts number of atoms to mass in grams. Supports both Exprs or SymPy quantities.

        Parameters
        ----------
        number : float or sympy.core.expr.Expr
            The number of atoms of the nuclide to be converted.
        atomic_mass : float or sympy.core.expr.Expr
            Atomic mass of the nuclide in g/mol.

        Returns
        -------
        float or sympy.core.expr.Expr
            Mass of material in grams.

        """

        return number / cls.avogadro * atomic_mass

    @classmethod
    def number_to_moles(cls, number: Union[float, Expr]) -> Union[float, Expr]:
        """
        Converts number of atoms to moles of nuclide.

        Parameters
        ----------
        number : float or sympy.core.expr.Expr
            The number of atoms of the nuclide to be converted.

        Returns
        -------
        float or sympy.core.expr.Expr
            Moles of nuclide.

        """

        return number / cls.avogadro


class QuantityConverterFloat(QuantityConverter):
    """
    Quantity converter using SymPy arbitrary precision operations.
    """

    avogadro: float = AVOGADRO

    @classmethod
    def __repr__(cls) -> str:
        return f"{cls.__name__} using double-precision floats."


class QuantityConverterSympy(QuantityConverter):
    """
    Quantity converter using SymPy arbitrary precision operations.
    """

    avogadro: Expr = nsimplify(AVOGADRO)

    @classmethod
    def __repr__(cls) -> str:
        return f"{cls.__name__} using SymPy arbitrary precision calculations."
