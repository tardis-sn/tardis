"""
Unit tests for converts.py classes and methods.
"""

import unittest

from sympy import Integer, log

from radioactivedecay.converters import (
    AVOGADRO,
    QuantityConverterFloat,
    QuantityConverterSympy,
    UnitConverterFloat,
    UnitConverterSympy,
)


class TestConverters(unittest.TestCase):
    """
    Unit tests for the converters.py constants.
    """

    def test_avogadro(self) -> None:
        """
        Test module constants.
        """

        self.assertEqual(AVOGADRO, 6.02214076e23)


class TestUnitConverterFloat(unittest.TestCase):
    """
    Unit tests for the converters.py UnitConverterFloat class.
    """

    def test_class_attributes(self) -> None:
        """
        Test class attributes of UnitConverterFloat objects.
        """

        self.assertEqual(UnitConverterFloat.time_units["s"], 1.0)
        self.assertEqual(UnitConverterFloat.time_units["y"], 86400.0)
        self.assertEqual(UnitConverterFloat.activity_units["Bq"], 1.0)
        self.assertEqual(UnitConverterFloat.mass_units["g"], 1.0)
        self.assertEqual(UnitConverterFloat.moles_units["mol"], 1.0)

        self.assertEqual(
            UnitConverterFloat.time_unit_err_msg,
            'is not a valid time unit, e.g. "s", "m", "h", "d" or "y".',
        )
        self.assertEqual(
            UnitConverterFloat.activity_unit_err_msg,
            'is not a valid activitiy unit, e.g. "Bq", "kBq", "Ci"...',
        )
        self.assertEqual(
            UnitConverterFloat.mass_unit_err_msg,
            'is not a valid mass unit, e.g. "g", "kg", "mg"...',
        )
        self.assertEqual(
            UnitConverterFloat.moles_unit_err_msg,
            'is not a valid moles unit, e.g. "mol", "kmol", "mmol"...',
        )

    def test_time_unit_conv_seconds(self) -> None:
        """
        Test conversion between seconds and different time units.
        """

        year_conv = 365.2422

        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "s", "s", year_conv), 1.0e0
        )
        self.assertAlmostEqual(
            UnitConverterFloat.time_unit_conv(1.0, "s", "ps", year_conv),
            1.0e12,
            places=(15 - 12),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.time_unit_conv(1.0, "s", "ns", year_conv),
            1.0e9,
            places=(15 - 9),
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "s", "μs", year_conv), 1.0e6
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "s", "us", year_conv), 1.0e6
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "s", "ms", year_conv), 1.0e3
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "s", "m", year_conv), 1.0 / 60.0
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "s", "h", year_conv),
            1.0 / (60.0**2),
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "s", "d", year_conv),
            1.0 / (60.0**2 * 24.0),
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "s", "y", year_conv),
            1.0 / (60.0**2 * 24.0 * year_conv),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.time_unit_conv(1.0, "ps", "s", year_conv),
            1.0e-12,
            places=(12 + 15),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.time_unit_conv(1.0, "ns", "s", year_conv),
            1.0e-9,
            places=(9 + 15),
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "μs", "s", year_conv), 1.0e-6
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "us", "s", year_conv), 1.0e-6
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "ms", "s", year_conv), 1.0e-3
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "m", "s", year_conv), 60.0
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "h", "s", year_conv), (60.0**2)
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "d", "s", year_conv),
            (60.0**2 * 24.0),
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "y", "s", year_conv),
            (60.0**2 * 24.0 * year_conv),
        )

        # Catch some incorrect time units
        with self.assertRaises(ValueError):
            UnitConverterFloat.time_unit_conv(1.0, "ty", "y", year_conv)
        with self.assertRaises(ValueError):
            UnitConverterFloat.time_unit_conv(1.0, "y", "ty", year_conv)
        with self.assertRaises(ValueError):
            UnitConverterFloat.time_unit_conv(1.0, "ty", 1.0, year_conv)

    def test_time_unit_conv_spelling_variations(self) -> None:
        """
        Test spelling variations of different time units.
        """

        year_conv = 365.2422

        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "s", "sec", year_conv), 1.0e0
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "s", "second", year_conv), 1.0e0
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "s", "seconds", year_conv), 1.0e0
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "h", "hr", year_conv), 1.0e0
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "h", "hour", year_conv), 1.0e0
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "h", "hours", year_conv), 1.0e0
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "d", "day", year_conv), 1.0e0
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "d", "days", year_conv), 1.0e0
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "y", "yr", year_conv), 1.0e0
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "y", "year", year_conv), 1.0e0
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "y", "years", year_conv), 1.0e0
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "sec", "s", year_conv), 1.0e0
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "second", "s", year_conv), 1.0e0
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "seconds", "s", year_conv), 1.0e0
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "hr", "h", year_conv), 1.0e0
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "hour", "h", year_conv), 1.0e0
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "hours", "h", year_conv), 1.0e0
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "day", "d", year_conv), 1.0e0
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "days", "d", year_conv), 1.0e0
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "yr", "y", year_conv), 1.0e0
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "year", "y", year_conv), 1.0e0
        )
        self.assertEqual(
            UnitConverterFloat.time_unit_conv(1.0, "years", "y", year_conv), 1.0e0
        )

    def test_time_unit_conv_year_prefixes(self) -> None:
        """
        Test conversions between different year prefixes.
        """

        year_conv = 365.2422

        self.assertAlmostEqual(
            UnitConverterFloat.time_unit_conv(1.0, "y", "ky", year_conv),
            1.0e-3,
            places=(3 + 15),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.time_unit_conv(1.0, "y", "My", year_conv),
            1.0e-6,
            places=(6 + 15),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.time_unit_conv(1.0, "y", "By", year_conv),
            1.0e-9,
            places=(9 + 15),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.time_unit_conv(1.0, "y", "Gy", year_conv),
            1.0e-9,
            places=(9 + 15),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.time_unit_conv(1.0, "y", "Ty", year_conv),
            1.0e-12,
            places=(12 + 15),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.time_unit_conv(1.0, "y", "Py", year_conv),
            1.0e-15,
            places=(15 + 15),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.time_unit_conv(1.0, "ky", "y", year_conv),
            1.0e3,
            places=(15 - 3),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.time_unit_conv(1.0, "My", "y", year_conv),
            1.0e6,
            places=(15 - 6),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.time_unit_conv(1.0, "Gy", "y", year_conv),
            1.0e9,
            places=(15 - 9),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.time_unit_conv(1.0, "Ty", "y", year_conv),
            1.0e12,
            places=(15 - 12),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.time_unit_conv(1.0, "Py", "y", year_conv),
            1.0e15,
            places=(15 - 15),
        )

    def test_activity_unit_conv(self) -> None:
        """
        Test conversions between activity units.
        """

        self.assertEqual(UnitConverterFloat.activity_unit_conv(1.0, "Bq", "Bq"), 1.0e0)
        self.assertAlmostEqual(
            UnitConverterFloat.activity_unit_conv(1.0, "Bq", "pBq"),
            1.0e12,
            places=(15 - 12),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.activity_unit_conv(1.0, "Bq", "nBq"),
            1.0e9,
            places=(15 - 9),
        )
        self.assertEqual(UnitConverterFloat.activity_unit_conv(1.0, "Bq", "μBq"), 1.0e6)
        self.assertEqual(UnitConverterFloat.activity_unit_conv(1.0, "Bq", "uBq"), 1.0e6)
        self.assertEqual(UnitConverterFloat.activity_unit_conv(1.0, "Bq", "mBq"), 1.0e3)
        self.assertAlmostEqual(
            UnitConverterFloat.activity_unit_conv(1.0, "Bq", "kBq"),
            1.0e-3,
            places=(3 + 15),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.activity_unit_conv(1.0, "Bq", "MBq"),
            1.0e-6,
            places=(6 + 15),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.activity_unit_conv(1.0, "Bq", "GBq"),
            1.0e-9,
            places=(9 + 15),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.activity_unit_conv(1.0, "Bq", "TBq"),
            1.0e-12,
            places=(12 + 15),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.activity_unit_conv(1.0, "Bq", "PBq"),
            1.0e-15,
            places=(15 + 15),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.activity_unit_conv(1.0, "Bq", "EBq"),
            1.0e-18,
            places=(18 + 15),
        )

        self.assertAlmostEqual(
            UnitConverterFloat.activity_unit_conv(1.0, "Bq", "pCi"),
            1.0e12 / 3.7e10,
            places=(15 - 12),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.activity_unit_conv(1.0, "Bq", "nCi"),
            1.0e9 / 3.7e10,
            places=(15 - 9),
        )
        self.assertEqual(
            UnitConverterFloat.activity_unit_conv(1.0, "Bq", "μCi"), 1.0e6 / 3.7e10
        )
        self.assertEqual(
            UnitConverterFloat.activity_unit_conv(1.0, "Bq", "uCi"), 1.0e6 / 3.7e10
        )
        self.assertEqual(
            UnitConverterFloat.activity_unit_conv(1.0, "Bq", "mCi"), 1.0e3 / 3.7e10
        )
        self.assertEqual(
            UnitConverterFloat.activity_unit_conv(1.0, "Bq", "Ci"), 1.0 / 3.7e10
        )
        self.assertAlmostEqual(
            UnitConverterFloat.activity_unit_conv(1.0, "Bq", "kCi"),
            1.0e-3 / 3.7e10,
            places=(3 + 15),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.activity_unit_conv(1.0, "Bq", "MCi"),
            1.0e-6 / 3.7e10,
            places=(6 + 15),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.activity_unit_conv(1.0, "Bq", "GCi"),
            1.0e-9 / 3.7e10,
            places=(9 + 15),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.activity_unit_conv(1.0, "Bq", "TCi"),
            1.0e-12 / 3.7e10,
            places=(12 + 15),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.activity_unit_conv(1.0, "Bq", "PCi"),
            1.0e-15 / 3.7e10,
            places=(15 + 15),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.activity_unit_conv(1.0, "Bq", "ECi"),
            1.0e-18 / 3.7e10,
            places=(18 + 15),
        )

        self.assertEqual(UnitConverterFloat.activity_unit_conv(1.0, "Bq", "dpm"), 60.0)

        # Catch some incorrect activity units
        with self.assertRaises(ValueError):
            UnitConverterFloat.activity_unit_conv(1.0, "tBq", "Bq")
        with self.assertRaises(ValueError):
            UnitConverterFloat.activity_unit_conv(1.0, "Bq", "tBq")
        with self.assertRaises(ValueError):
            UnitConverterFloat.activity_unit_conv(1.0, "tBq", 1.0)

    def test_mass_unit_conv(self) -> None:
        """
        Test conversions between mass units.
        """

        self.assertEqual(UnitConverterFloat.mass_unit_conv(1.0, "g", "g"), 1.0e0)
        self.assertAlmostEqual(
            UnitConverterFloat.mass_unit_conv(1.0, "g", "pg"), 1.0e12, places=(15 - 12)
        )
        self.assertAlmostEqual(
            UnitConverterFloat.mass_unit_conv(1.0, "g", "ng"), 1.0e9, places=(15 - 9)
        )
        self.assertEqual(UnitConverterFloat.mass_unit_conv(1.0, "g", "μg"), 1.0e6)
        self.assertEqual(UnitConverterFloat.mass_unit_conv(1.0, "g", "ug"), 1.0e6)
        self.assertEqual(UnitConverterFloat.mass_unit_conv(1.0, "g", "mg"), 1.0e3)
        self.assertAlmostEqual(
            UnitConverterFloat.mass_unit_conv(1.0, "g", "kg"), 1.0e-3, places=(3 + 15)
        )
        self.assertAlmostEqual(
            UnitConverterFloat.mass_unit_conv(1.0, "g", "Mg"), 1.0e-6, places=(6 + 15)
        )
        self.assertAlmostEqual(
            UnitConverterFloat.mass_unit_conv(1.0, "g", "t"), 1.0e-6, places=(6 + 15)
        )
        self.assertAlmostEqual(
            UnitConverterFloat.mass_unit_conv(1.0, "g", "ton"), 1.0e-6, places=(6 + 15)
        )

        # Catch some incorrect activity units
        with self.assertRaises(ValueError):
            UnitConverterFloat.mass_unit_conv(1.0, "tg", "g")
        with self.assertRaises(ValueError):
            UnitConverterFloat.mass_unit_conv(1.0, "g", "tg")
        with self.assertRaises(ValueError):
            UnitConverterFloat.mass_unit_conv(1.0, "tg", 1.0)

    def test_moles_unit_conv(self) -> None:
        """
        Test conversions between moles orders of magnitude.
        """

        self.assertEqual(UnitConverterFloat.moles_unit_conv(1.0, "mol", "mol"), 1.0e0)
        self.assertAlmostEqual(
            UnitConverterFloat.moles_unit_conv(1.0, "mol", "pmol"),
            1.0e12,
            places=(15 - 12),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.moles_unit_conv(1.0, "mol", "nmol"),
            1.0e9,
            places=(15 - 9),
        )
        self.assertEqual(UnitConverterFloat.moles_unit_conv(1.0, "mol", "μmol"), 1.0e6)
        self.assertEqual(UnitConverterFloat.moles_unit_conv(1.0, "mol", "umol"), 1.0e6)
        self.assertEqual(UnitConverterFloat.moles_unit_conv(1.0, "mol", "mmol"), 1.0e3)
        self.assertAlmostEqual(
            UnitConverterFloat.moles_unit_conv(1.0, "mol", "kmol"),
            1.0e-3,
            places=(3 + 15),
        )
        self.assertAlmostEqual(
            UnitConverterFloat.moles_unit_conv(1.0, "mol", "Mmol"),
            1.0e-6,
            places=(6 + 15),
        )

        # Catch some incorrect activity units
        with self.assertRaises(ValueError):
            UnitConverterFloat.moles_unit_conv(1.0, "tmol", "mol")
        with self.assertRaises(ValueError):
            UnitConverterFloat.moles_unit_conv(1.0, "mol", "tmol")
        with self.assertRaises(ValueError):
            UnitConverterFloat.moles_unit_conv(1.0, "tmol", 1.0)

    def test___repr__(self) -> None:
        """
        Test UnitConverterFloat __repr__ strings.
        """

        self.assertEqual(
            repr(UnitConverterFloat()),
            "UnitConverterFloat using double-precision floats.",
        )


class TestUnitConverterSympy(unittest.TestCase):
    """
    Unit tests for the converters.py UnitConverterSympy class.
    """

    def test_class_attributes(self) -> None:
        """
        Test class attributes of UnitConverterSympy objects.
        """

        self.assertEqual(UnitConverterSympy.time_units["s"], Integer(1))
        self.assertEqual(UnitConverterSympy.time_units["y"], Integer(86400))
        self.assertEqual(UnitConverterSympy.activity_units["Bq"], Integer(1))
        self.assertEqual(UnitConverterSympy.mass_units["g"], Integer(1))
        self.assertEqual(UnitConverterSympy.moles_units["mol"], Integer(1))

    def test_time_unit_conv(self) -> None:
        """
        Test of the SymPy version of time_unit_conv().
        """

        year_conv = Integer(3652422) / 10000

        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "ps", "ns", year_conv),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "ns", "us", year_conv),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "μs", "ms", year_conv),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "us", "ms", year_conv),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "ms", "s", year_conv),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "s", "m", year_conv),
            1 / Integer(60),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "m", "h", year_conv),
            1 / Integer(60),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "h", "d", year_conv),
            1 / Integer(24),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "d", "y", year_conv),
            10000 / Integer(3652422),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "y", "ky", year_conv),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "ky", "My", year_conv),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "My", "By", year_conv),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "My", "Gy", year_conv),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "By", "Ty", year_conv),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "Gy", "Ty", year_conv),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "Ty", "Py", year_conv),
            1 / Integer(1000),
        )

        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "ns", "ps", year_conv),
            Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "μs", "ns", year_conv),
            Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "us", "ns", year_conv),
            Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "ms", "us", year_conv),
            Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "s", "ms", year_conv),
            Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "m", "s", year_conv),
            Integer(60),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "h", "m", year_conv),
            Integer(60),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "d", "h", year_conv),
            Integer(24),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "y", "d", year_conv),
            Integer(3652422) / 10000,
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "ky", "y", year_conv),
            Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "My", "ky", year_conv),
            Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "By", "My", year_conv),
            Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "Gy", "My", year_conv),
            Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "Ty", "By", year_conv),
            Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "Ty", "Gy", year_conv),
            Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.time_unit_conv(Integer(1), "Py", "Ty", year_conv),
            Integer(1000),
        )

        # Catch some incorrect time units
        with self.assertRaises(ValueError):
            UnitConverterSympy.time_unit_conv(Integer(1), "ty", "y", year_conv)
        with self.assertRaises(ValueError):
            UnitConverterSympy.time_unit_conv(Integer(1), "y", "ty", year_conv)
        with self.assertRaises(ValueError):
            UnitConverterSympy.time_unit_conv(Integer(1), "ty", Integer(1), year_conv)

    def test_activity_unit_conv(self) -> None:
        """
        Test of the SymPy version of activity_unit_conv().
        """

        self.assertEqual(
            UnitConverterSympy.activity_unit_conv(Integer(1), "pBq", "nBq"),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.activity_unit_conv(Integer(1), "nBq", "μBq"),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.activity_unit_conv(Integer(1), "μBq", "uBq"), Integer(1)
        )
        self.assertEqual(
            UnitConverterSympy.activity_unit_conv(Integer(1), "uBq", "mBq"),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.activity_unit_conv(Integer(1), "mBq", "Bq"),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.activity_unit_conv(Integer(1), "Bq", "kBq"),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.activity_unit_conv(Integer(1), "kBq", "MBq"),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.activity_unit_conv(Integer(1), "MBq", "GBq"),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.activity_unit_conv(Integer(1), "GBq", "TBq"),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.activity_unit_conv(Integer(1), "TBq", "PBq"),
            1 / Integer(1000),
        )

        self.assertEqual(
            UnitConverterSympy.activity_unit_conv(Integer(1), "pBq", "pCi"),
            1 / Integer(37000000000),
        )
        self.assertEqual(
            UnitConverterSympy.activity_unit_conv(Integer(1), "pCi", "nCi"),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.activity_unit_conv(Integer(1), "nCi", "μCi"),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.activity_unit_conv(Integer(1), "μCi", "uCi"), Integer(1)
        )
        self.assertEqual(
            UnitConverterSympy.activity_unit_conv(Integer(1), "uCi", "mCi"),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.activity_unit_conv(Integer(1), "mCi", "Ci"),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.activity_unit_conv(Integer(1), "Ci", "kCi"),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.activity_unit_conv(Integer(1), "kCi", "MCi"),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.activity_unit_conv(Integer(1), "MCi", "GCi"),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.activity_unit_conv(Integer(1), "GCi", "TCi"),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.activity_unit_conv(Integer(1), "TCi", "PCi"),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.activity_unit_conv(Integer(1), "PCi", "ECi"),
            1 / Integer(1000),
        )

        self.assertEqual(
            UnitConverterSympy.activity_unit_conv(Integer(1), "Bq", "dpm"),
            Integer(60),
        )

        # Catch some incorrect activity units
        with self.assertRaises(ValueError):
            UnitConverterSympy.activity_unit_conv(Integer(1), "tBq", "Bq")
        with self.assertRaises(ValueError):
            UnitConverterSympy.activity_unit_conv(Integer(1), "Bq", "tBq")
        with self.assertRaises(ValueError):
            UnitConverterSympy.activity_unit_conv(Integer(1), "tBq", Integer(1))

    def test_mass_unit_conv(self) -> None:
        """
        Test of the SymPy version of mass_unit_conv().
        """

        self.assertEqual(
            UnitConverterSympy.mass_unit_conv(Integer(1), "pg", "ng"), 1 / Integer(1000)
        )
        self.assertEqual(
            UnitConverterSympy.mass_unit_conv(Integer(1), "ng", "μg"), 1 / Integer(1000)
        )
        self.assertEqual(
            UnitConverterSympy.mass_unit_conv(Integer(1), "μg", "ug"), Integer(1)
        )
        self.assertEqual(
            UnitConverterSympy.mass_unit_conv(Integer(1), "ug", "mg"), 1 / Integer(1000)
        )
        self.assertEqual(
            UnitConverterSympy.mass_unit_conv(Integer(1), "mg", "g"), 1 / Integer(1000)
        )
        self.assertEqual(
            UnitConverterSympy.mass_unit_conv(Integer(1), "g", "kg"), 1 / Integer(1000)
        )
        self.assertEqual(
            UnitConverterSympy.mass_unit_conv(Integer(1), "kg", "Mg"), 1 / Integer(1000)
        )
        self.assertEqual(
            UnitConverterSympy.mass_unit_conv(Integer(1), "Mg", "t"), Integer(1)
        )
        self.assertEqual(
            UnitConverterSympy.mass_unit_conv(Integer(1), "Mg", "ton"), Integer(1)
        )

        # Catch some incorrect activity units
        with self.assertRaises(ValueError):
            UnitConverterSympy.mass_unit_conv(Integer(1), "tg", "g")
        with self.assertRaises(ValueError):
            UnitConverterSympy.mass_unit_conv(Integer(1), "g", "tg")
        with self.assertRaises(ValueError):
            UnitConverterSympy.mass_unit_conv(Integer(1), "tg", Integer(1))

    def test_moles_unit_conv(self) -> None:
        """
        Test of the SymPy version of moles_unit_conv().
        """

        self.assertEqual(
            UnitConverterSympy.moles_unit_conv(Integer(1), "pmol", "nmol"),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.moles_unit_conv(Integer(1), "nmol", "μmol"),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.moles_unit_conv(Integer(1), "μmol", "umol"), Integer(1)
        )
        self.assertEqual(
            UnitConverterSympy.moles_unit_conv(Integer(1), "umol", "mmol"),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.moles_unit_conv(Integer(1), "mmol", "mol"),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.moles_unit_conv(Integer(1), "mol", "kmol"),
            1 / Integer(1000),
        )
        self.assertEqual(
            UnitConverterSympy.moles_unit_conv(Integer(1), "kmol", "Mmol"),
            1 / Integer(1000),
        )

        # Catch some incorrect activity units
        with self.assertRaises(ValueError):
            UnitConverterSympy.moles_unit_conv(Integer(1), "tmol", "mol")
        with self.assertRaises(ValueError):
            UnitConverterSympy.moles_unit_conv(Integer(1), "mol", "tmol")
        with self.assertRaises(ValueError):
            UnitConverterSympy.moles_unit_conv(Integer(1), "tmol", Integer(1))

    def test___repr__(self) -> None:
        """
        Test UnitConverterSympy __repr__ strings.
        """

        self.assertEqual(
            repr(UnitConverterSympy()),
            "UnitConverterSympy using SymPy arbitrary precision calculations.",
        )


class TestQuantityConverterFloat(unittest.TestCase):
    """
    Unit tests for the converters.py QuantityConverterFloat class.
    """

    def test_class_attribute(self) -> None:
        """
        Test class attribute of QuantityConverterFloat.
        """

        self.assertEqual(QuantityConverterFloat.avogadro, 6.02214076e23)

    def test_activity_to_number(self) -> None:
        """
        Test the conversion of activity in Bq to number of atoms.
        """

        self.assertEqual(
            QuantityConverterFloat.activity_to_number(1.0, 1.7828715741004621e-09),
            560892895.7794082,
        )

    def test_mass_to_number(self) -> None:
        """
        Test the conversion of mass in grams to number of atoms.
        """

        self.assertEqual(
            QuantityConverterFloat.mass_to_number(1.0, 3.01604928132),
            1.996698395247825e23,
        )
        self.assertEqual(
            QuantityConverterFloat.mass_to_number(1.0, 3.01602932197),
            1.9967116089131645e23,
        )

    def test_moles_to_number(self) -> None:
        """
        Test the conversion of number of moles to number of atoms.
        """

        self.assertEqual(QuantityConverterFloat.moles_to_number(1.0), 6.02214076e23)
        self.assertEqual(QuantityConverterFloat.moles_to_number(0.0), 0.0)

    def test_number_to_activity(self) -> None:
        """
        Test the conversion of number of atoms to activity in Bq.
        """

        self.assertEqual(
            QuantityConverterFloat.number_to_activity(
                560892895.7794082, 1.7828715741004621e-09
            ),
            1.0,
        )
        self.assertEqual(QuantityConverterFloat.number_to_activity(1.0, 0.0), 0.0)
        self.assertEqual(QuantityConverterFloat.number_to_activity(0.0, 0.0), 0.0)

    def test_number_to_mass(self) -> None:
        """
        Test the conversion of number of atoms to mass in grams.
        """

        self.assertEqual(
            QuantityConverterFloat.number_to_mass(1.996698395247825e23, 3.01604928132),
            1.0,
        )
        self.assertEqual(
            QuantityConverterFloat.number_to_mass(1.9967116089131645e23, 3.01602932197),
            1.0,
        )

    def test_number_to_moles(self) -> None:
        """
        Test the conversion of number of atoms to number of moles.
        """

        self.assertEqual(QuantityConverterFloat.number_to_moles(6.02214076e23), 1.0)
        self.assertEqual(QuantityConverterFloat.number_to_moles(0.0), 0.0)

    def test___repr__(self) -> None:
        """
        Test QuantityConverterFloat __repr__ strings.
        """

        self.assertEqual(
            repr(QuantityConverterFloat()),
            "QuantityConverterFloat using double-precision floats.",
        )


class TestQuantityConverterSympy(unittest.TestCase):
    """
    Unit tests for the converters.py QuantityConverterSympy class.
    """

    def test_class_attribute(self) -> None:
        """
        Test class attribute of QuantityConverterSympy class.
        """

        self.assertEqual(
            QuantityConverterSympy.avogadro, Integer(602214076000000000000000)
        )

    def test_activity_to_number(self) -> None:
        """
        Test the conversion of activity in Bq to number of atoms.
        """

        self.assertEqual(
            QuantityConverterSympy.activity_to_number(
                Integer(1), (Integer(625) * log(2)) / Integer(242988330816)
            ),
            Integer(242988330816) / (Integer(625) * log(2)),
        )

    def test_mass_to_number(self) -> None:
        """
        Test the conversion of mass in grams to number of atoms.
        """

        self.assertEqual(
            QuantityConverterSympy.mass_to_number(
                Integer(1),
                Integer("15055351900000000000000000000000000") / Integer(75401232033),
            ),
            Integer(75401232033)
            / Integer("15055351900000000000000000000000000")
            * Integer("602214076000000000000000"),
        )
        self.assertEqual(
            QuantityConverterSympy.mass_to_number(
                Integer(1),
                Integer("60221407600000000000000000000000000") / Integer(301602932197),
            ),
            Integer(301602932197)
            / Integer("60221407600000000000000000000000000")
            * Integer("602214076000000000000000"),
        )

    def test_moles_to_number(self) -> None:
        """
        Test the conversion of number of moles to number of atoms.
        """

        self.assertEqual(
            QuantityConverterSympy.moles_to_number(Integer(1)),
            Integer("602214076000000000000000"),
        )
        self.assertEqual(QuantityConverterSympy.moles_to_number(Integer(0)), Integer(0))

    def test_number_to_activity(self) -> None:
        """
        Test the conversion of number of atoms to activity in Bq.
        """

        self.assertEqual(
            QuantityConverterSympy.number_to_activity(
                Integer(1), (Integer(625) * log(2)) / Integer(242988330816)
            ),
            Integer(625) * log(2) / Integer(242988330816),
        )
        self.assertEqual(
            QuantityConverterSympy.number_to_activity(
                Integer(0), (Integer(625) * log(2)) / Integer(242988330816)
            ),
            Integer(0),
        )
        self.assertEqual(
            QuantityConverterSympy.number_to_activity(Integer(0), Integer(0)),
            Integer(0),
        )

    def test_number_to_mass(self) -> None:
        """
        Test the conversion of number of atoms to mass in grams.
        """

        self.assertEqual(
            QuantityConverterSympy.number_to_mass(
                Integer(1),
                Integer("15055351900000000000000000000000000") / Integer(75401232033),
            ),
            Integer("15055351900000000000000000000000000")
            / Integer(75401232033)
            / Integer("602214076000000000000000"),
        )
        self.assertEqual(
            QuantityConverterSympy.number_to_mass(
                Integer(1),
                Integer("60221407600000000000000000000000000") / Integer(301602932197),
            ),
            Integer("60221407600000000000000000000000000")
            / Integer(301602932197)
            / Integer("602214076000000000000000"),
        )

    def test_number_to_moles(self) -> None:
        """
        Test the conversion of number of atoms to number of moles.
        """

        self.assertEqual(
            QuantityConverterSympy.number_to_moles(Integer("602214076000000000000000")),
            Integer(1),
        )
        self.assertEqual(QuantityConverterSympy.number_to_moles(Integer(0)), Integer(0))

    def test___repr__(self) -> None:
        """
        Test QuantityConveterSympy __repr__ strings.
        """

        self.assertEqual(
            repr(QuantityConverterSympy()),
            "QuantityConverterSympy using SymPy arbitrary precision calculations.",
        )


if __name__ == "__main__":
    unittest.main()
