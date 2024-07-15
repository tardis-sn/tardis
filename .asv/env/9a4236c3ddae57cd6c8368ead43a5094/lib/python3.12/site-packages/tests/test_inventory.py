"""
Unit tests for inventory.py classes and methods.
"""

import copy
import math
import unittest
import warnings
from typing import Dict
from unittest.mock import patch

import numpy as np
import pandas as pd
from sympy import Integer, log

from radioactivedecay.decaydata import DEFAULTDATA, load_dataset
from radioactivedecay.inventory import Inventory, InventoryHP
from radioactivedecay.nuclide import Nuclide

# pylint: disable=protected-access, too-many-public-methods


np.set_printoptions(legacy="1.25")


def warning_message_if_dict_not_equal(
    calculated: Dict[str, float], expected: Dict[str, float]
) -> None:
    """
    Warning message if calculated dictionary of floats is not equal to expected dictionary of
    floats.
    """

    return (
        f"Warning: calculated result:\n{calculated}\nis not identical to expected result:\n"
        f"{expected}\nFloating point error build up or other code issue has occurred."
    )


def dict_assert_almost_equal(
    test_case: unittest.TestCase,
    dict_a: Dict[str, float],
    dict_b: Dict[str, float],
) -> None:
    """
    Check whether two dictionaries with str keys and float values have:
        i) the same keys
        ii) float values are close (using `math.isclose()`)

    This function could be improved by outputting a better traceback in the event that an
    assertTrue fails.
    """

    test_case.assertEqual(set(dict_a), set(dict_b))

    for key in dict_a:
        try:
            test_case.assertTrue(
                math.isclose(dict_a[key], dict_b[key], rel_tol=1e-7, abs_tol=1e-30)
            )
        except AssertionError as error:
            raise AssertionError(
                f"Floats are not within error tolerance: {dict_a[key]} {dict_b[key]} for nuclide "
                f"{key}."
            ) from error


class TestInventory(unittest.TestCase):
    """
    Unit tests for the inventory.py Inventory class.
    """

    def test_instantiation(self) -> None:
        """
        Test instantiation of Inventory objects.
        """

        # check basic instantiation
        inv = Inventory({"H3": 1}, "num", True)
        self.assertEqual(inv.contents, {"H-3": 1})
        self.assertEqual(inv.decay_data, DEFAULTDATA)
        self.assertEqual(inv.decay_matrices, DEFAULTDATA.scipy_data)

        # check instantiation using Nuclide instance
        tritium = Nuclide("H3")
        inv = Inventory({tritium: 1.0}, "num", True)
        self.assertEqual(inv.contents, {"H-3": 1})

        # check that check=False turns off nuclide string conversion and quantity check
        inv = Inventory({"H3": 1.0}, "num", False)
        self.assertEqual(inv.contents, {"H3": 1.0})
        inv = Inventory({"H-3": -1.0}, "num", False)
        self.assertEqual(inv.contents, {"H-3": -1.0})

        # check sorting of contents dictionary
        inv = Inventory({"He-3": 2.0, "H-3": 1.0}, "num", False)
        self.assertEqual(inv.contents, {"H-3": 1.0, "He-3": 2.0})

        # check conversion of input quantity to number
        inv = Inventory({"H-3": 1.0}, "mBq")
        self.assertEqual(inv.contents, {"H-3": 560892.8957794083})

        # check instantiation with a stable nuclide activity raises a ValueError
        with self.assertRaises(ValueError):
            Inventory({"He-3": 0.0}, "Bq", False)

        # check instantiation with NumPy int datatype
        inv = Inventory({"H-3": np.int32(1)}, "num")
        self.assertEqual(inv.contents, {"H-3": 1})

    def test__parse_nuclides(self) -> None:
        """
        Test the conversion of nuclide strings or Nuclide instances in a contents dictionary
        to nuclide strings into Ab-XY form.
        """

        nuclides = DEFAULTDATA.nuclides
        dataset_name = DEFAULTDATA.dataset_name
        self.assertEqual(
            Inventory._parse_nuclides({"H-3": 1.0}, nuclides, dataset_name),
            {"H-3": 1.0},
        )
        self.assertEqual(
            Inventory._parse_nuclides({"He3": 1.0}, nuclides, dataset_name),
            {"He-3": 1.0},
        )

        tritium = Nuclide("H-3")
        self.assertEqual(
            Inventory._parse_nuclides(
                {tritium: 1.0, "3He": 2.0}, nuclides, dataset_name
            ),
            {"H-3": 1.0, "He-3": 2.0},
        )

    def test__check_values(self) -> None:
        """
        Test that the contents dictionary values (amounts of each nuclide) are physical.
        """

        self.assertIsNone(Inventory._check_values({"H-3": 1.0}))
        self.assertIsNone(Inventory._check_values({"H-3": 0}))
        self.assertIsNone(Inventory._check_values({"H-3": Integer(0)}))

        with self.assertRaises(ValueError):
            Inventory._check_values({"H-3": -1})
        with self.assertRaises(ValueError):
            Inventory._check_values({"H-3": "1"})

    def test__convert_to_number(self) -> None:
        """
        Test that the contents dictionary values (amounts of each nuclide) are physical.
        """

        inv = Inventory({"H3": 1}, "num", True)
        self.assertEqual(
            inv._convert_to_number({"H-3": 1.0}, "num", "dataset_name"),
            {"H-3": 1.0},
        )
        self.assertEqual(
            inv._convert_to_number({"H-3": 1.0}, "Bq", "dataset_name"),
            {"H-3": 560892895.7794082},
        )
        self.assertEqual(
            inv._convert_to_number({"H-3": 1.0}, "mBq", "dataset_name"),
            {"H-3": 560892.8957794083},
        )
        self.assertEqual(
            inv._convert_to_number({"H-3": 1.0}, "mol", "dataset_name"),
            {"H-3": 6.02214076e23},
        )
        self.assertEqual(
            inv._convert_to_number({"He-3": 1.0}, "kg", "dataset_name"),
            {"He-3": 1.9967116089131648e26},
        )

        # check catch of incorrect units
        with self.assertRaises(ValueError):
            inv._convert_to_number({"H-3": 1.0}, "xyz", "dataset_name")

    def test_nuclides(self) -> None:
        """
        Test Inventory.nuclides property.
        """

        inv = Inventory({"H-3": 1.0})
        self.assertEqual(inv.nuclides, ["H-3"])
        inv = Inventory({"Tc-99m": 2.3, "I-123": 5.8})
        self.assertEqual(inv.nuclides, ["I-123", "Tc-99m"])

    def test_numbers(self) -> None:
        """
        Test Inventory.numbers() method.
        """

        inv = Inventory({"H-3": 1}, "num")
        self.assertEqual(inv.numbers(), {"H-3": 1})
        inv = Inventory({"Tc-99m": 2.3, "I-123": 5.8})
        self.assertEqual(
            inv.numbers(), {"I-123": 399738.47946141585, "Tc-99m": 71852.27235544211}
        )

    def test_activities(self) -> None:
        """
        Test Inventory.activities() method.
        """

        inv = Inventory({"H-3": 1})
        self.assertEqual(inv.activities(), {"H-3": 1})
        inv = Inventory({"Tc-99m": 2.3, "I-123": 5.8})
        self.assertEqual(inv.activities(), {"I-123": 5.8, "Tc-99m": 2.3})
        inv = Inventory({"Tc-99m": 2300, "I-123": 5800})
        self.assertEqual(inv.activities("kBq"), {"I-123": 5.8, "Tc-99m": 2.3})

    def test_masses(self) -> None:
        """
        Test Inventory.masses() method.
        """

        inv = Inventory({"H-3": 1}, "g")
        self.assertEqual(inv.masses(), {"H-3": 1})
        inv = Inventory({"Tc-99m": 2.3, "I-123": 5.8})
        self.assertEqual(
            inv.masses(),
            {"I-123": 8.158243973887584e-17, "Tc-99m": 1.1800869622748501e-17},
        )
        self.assertEqual(
            inv.masses("pg"),
            {"I-123": 8.158243973887584e-05, "Tc-99m": 1.1800869622748502e-05},
        )

    def test_moles(self) -> None:
        """
        Test Inventory.moles() method.
        """

        inv = Inventory({"H-3": 1}, "mol")
        self.assertEqual(inv.moles(), {"H-3": 1})
        inv = Inventory({"Tc-99m": 2.3, "I-123": 5.8})
        self.assertEqual(
            inv.moles(),
            {"I-123": 6.637813617983513e-19, "Tc-99m": 1.1931350531142702e-19},
        )
        self.assertEqual(
            inv.moles("mmol"),
            {"I-123": 6.637813617983513e-16, "Tc-99m": 1.1931350531142702e-16},
        )

    def test_activity_fractions(self) -> None:
        """
        Test Inventory.activity_fractions() method.
        """

        inv = Inventory({"Tc-99m": 2.3, "I-123": 5.8})
        self.assertEqual(
            inv.activity_fractions(),
            {"I-123": 0.7160493827160493, "Tc-99m": 0.2839506172839506},
        )
        self.assertEqual(sum(inv.activity_fractions().values()), 1.0)

    def test_mass_fractions(self) -> None:
        """
        Test Inventory.mass_fractions() method.
        """

        inv = Inventory({"Tc-99m": 2.3, "I-123": 5.8})
        self.assertEqual(
            inv.mass_fractions(),
            {"I-123": 0.8736297770616593, "Tc-99m": 0.12637022293834066},
        )
        self.assertEqual(sum(inv.mass_fractions().values()), 1.0)

    def test_mole_fractions(self) -> None:
        """
        Test Inventory.mole_fractions() method.
        """

        inv = Inventory({"Tc-99m": 2.3, "I-123": 5.8})
        self.assertEqual(
            inv.mole_fractions(),
            {"I-123": 0.8476385041932588, "Tc-99m": 0.15236149580674116},
        )
        self.assertEqual(sum(inv.mole_fractions().values()), 1.0)

    def test___len__(self) -> None:
        """
        Test len() on Inventory.
        """

        inv = Inventory({"H-3": 1})
        self.assertEqual(len(inv), 1)
        inv = Inventory({"Tc-99m": 2.3, "I-123": 5.8})
        self.assertEqual(len(inv), 2)

    def test_add(self) -> None:
        """
        Test Inventory.add() method to append to an inventory.
        """

        inv = Inventory({"H-3": 1}, "num")
        inv.add({60140000: 3.0, "K-40": 4.0}, "num")
        self.assertEqual(inv.contents, {"C-14": 3.0, "H-3": 1.0, "K-40": 4.0})
        inv.add({"H-3": 3.0}, "num")
        self.assertEqual(inv.contents, {"C-14": 3.0, "H-3": 4.0, "K-40": 4.0})

        inv = Inventory({"H-3": 1}, "num")
        inv.add({Nuclide("C-14"): 3.0, "K-40": 4.0}, "num")
        self.assertEqual(inv.contents, {"C-14": 3.0, "H-3": 1.0, "K-40": 4.0})
        inv.add({Nuclide("H-3"): 3.0}, "num")
        self.assertEqual(inv.contents, {"C-14": 3.0, "H-3": 4.0, "K-40": 4.0})

    def test_subtract(self) -> None:
        """
        Test Inventory.subtract() method to take away a dictionary from an inventory.
        """

        inv = Inventory({"C-14": 3.0, "H-3": 4.0, "K-40": 4.0}, "num")
        inv.subtract({60140000: 3.0, "K-40": 4.0}, "num")
        self.assertEqual(inv.contents, {"C-14": 0.0, "H-3": 4.0, "K-40": 0.0}, "num")

        inv = Inventory({"C-14": 3.0, "H-3": 4.0, "K-40": 4.0}, "num")
        inv.subtract({"C-14": 3.0, Nuclide("K-40"): 4.0}, "num")
        self.assertEqual(inv.contents, {"C-14": 0.0, "H-3": 4.0, "K-40": 0.0}, "num")

    def test___add__(self) -> None:
        """
        Test operator to add two inventory objects together.
        """

        inv1 = Inventory({"H-3": 1.0}, "num")
        inv2 = Inventory({"C-14": 1.0, "H-3": 4.0}, "num")
        inv = inv1 + inv2
        self.assertEqual(inv.contents, {"C-14": 1.0, "H-3": 5.0}, "num")

        temp_data = copy.deepcopy(DEFAULTDATA)
        temp_data.dataset_name = "icrp107_"
        inv3 = Inventory({"H-3": 2.0}, "num", decay_data=temp_data)
        with self.assertRaises(ValueError):
            inv = inv1 + inv3

    def test___subtract__(self) -> None:
        """
        Test operator to subtract one inventory object from another.
        """

        inv1 = Inventory({"H-3": 1.0}, "num")
        inv2 = Inventory({"C-14": 1.0, "H-3": 4.0}, "num")
        inv = inv2 - inv1
        self.assertEqual(inv.contents, {"C-14": 1.0, "H-3": 3.0})

        temp_data = copy.deepcopy(DEFAULTDATA)
        temp_data.dataset_name = "icrp107_"
        inv3 = Inventory({"H-3": 2.0}, "num", decay_data=temp_data)
        with self.assertRaises(ValueError):
            inv = inv1 - inv3

    def test___mul__(self) -> None:
        """
        Test operator to multiply activities in inventory by constant.
        """

        inv = Inventory({"Sr-90": 1.0, "Cs-137": 1.0}, "num")
        inv = inv * 2
        self.assertEqual(inv.contents, {"Cs-137": 2.0, "Sr-90": 2.0})

    def test___rmul__(self) -> None:
        """
        Test operator to right multiply constant by activities in inventory.
        """

        inv = Inventory({"Sr-90": 1.0, "Cs-137": 1.0}, "num")
        inv = 2 * inv
        self.assertEqual(inv.contents, {"Cs-137": 2.0, "Sr-90": 2.0})

    def test___truediv__(self) -> None:
        """
        Test operator to multiply activities in inventory by constant.
        """

        inv = Inventory({"Sr-90": 1.0, "Cs-137": 1.0}, "num")
        inv = inv / 2
        self.assertEqual(inv.contents, {"Cs-137": 0.5, "Sr-90": 0.5})

    def test_remove(self) -> None:
        """
        Test operator to remove nuclides from an inventory.
        """

        inv = Inventory({"C-14": 3.0, "H-3": 4.0, "K-40": 4.0}, "num")
        with self.assertRaises(NotImplementedError):
            inv.remove(1.0)

    def test_remove_string(self) -> None:
        """
        Test operator to remove one nuclide from an inventory using a nuclide string.
        """

        inv = Inventory({"C-14": 3.0, "H-3": 4.0, "K-40": 4.0}, "num")
        inv.remove("H-3")
        self.assertEqual(inv.contents, {"C-14": 3.0, "K-40": 4.0})

        with self.assertRaises(ValueError):
            inv.remove("Be-10")

    def test_remove_canonical_id(self) -> None:
        """
        Test operator to remove one nuclide from an inventory using a nuclide canonical id.
        """

        inv = Inventory({"C-14": 3.0, "H-3": 4.0, "K-40": 4.0}, "num")
        inv.remove(10030000)
        self.assertEqual(inv.contents, {"C-14": 3.0, "K-40": 4.0})

        with self.assertRaises(ValueError):
            inv.remove(40100000)

    def test_remove_nuclide(self) -> None:
        """
        Test operator to remove one nuclide from an inventory using a ``Nuclide`` object.
        """

        inv = Inventory({"C-14": 3.0, "H-3": 4.0, "K-40": 4.0}, "num")
        inv.remove(Nuclide("H-3"))
        self.assertEqual(inv.contents, {"C-14": 3.0, "K-40": 4.0})

        with self.assertRaises(ValueError):
            inv.remove(Nuclide("Be-10"))

    def test_remove_list(self) -> None:
        """
        Test operator to remove list of nuclides from an inventory.
        """

        inv = Inventory({"C-14": 3.0, "H-3": 4.0, "K-40": 4.0}, "num")
        inv.remove([10030000, "C-14"])
        self.assertEqual(inv.contents, {"K-40": 4.0})

        with self.assertRaises(ValueError):
            inv.remove(["Be-10", "C-14"])

        inv = Inventory({"C-14": 3.0, "H-3": 4.0, "K-40": 4.0}, "num")
        inv.remove(["H-3", Nuclide("C-14")])
        self.assertEqual(inv.contents, {"K-40": 4.0})

    def test_decay(self) -> None:
        """
        Test Inventory.decay() calculations.
        """

        inv = Inventory({"H-3": 10.0}, "Bq")
        self.assertEqual(inv.decay(12.32, "y").activities(), {"H-3": 5.0, "He-3": 0.0})

        inv = Inventory({"Tc-99m": 2.3, "I-123": 5.8}, "Bq")
        calculated = inv.decay(20.0, "h").activities()
        expected = {
            "I-123": 2.040459244534774,
            "Ru-99": 0.0,
            "Sb-123": 0.0,
            "Tc-99": 6.729944738772211e-09,
            "Tc-99m": 0.22950748010063513,
            "Te-123": 9.485166535243877e-18,
            "Te-123m": 7.721174031572363e-07,
        }
        if calculated != expected:
            warnings.warn(warning_message_if_dict_not_equal(calculated, expected))
        dict_assert_almost_equal(self, calculated, expected)

        inv = Inventory({"U-238": 99.274, "U-235": 0.720, "U-234": 0.005}, "Bq")
        calculated = inv.decay(1e9, "y").activities()
        expected = {
            "Ac-227": 0.2690006281740556,
            "At-218": 0.017002868638497183,
            "At-219": 2.227325201281319e-07,
            "Bi-210": 85.01434361515662,
            "Bi-211": 0.26900084425585846,
            "Bi-214": 85.01432618961896,
            "Bi-215": 2.1605054452429237e-07,
            "Fr-223": 0.0037122086688021884,
            "Hg-206": 1.6152725286830197e-06,
            "Pa-231": 0.2690006198549055,
            "Pa-234": 0.13601313171698984,
            "Pa-234m": 85.00820732310412,
            "Pb-206": 0.0,
            "Pb-207": 0.0,
            "Pb-210": 85.01434361489548,
            "Pb-211": 0.2690008442558569,
            "Pb-214": 84.99734032384839,
            "Po-210": 85.01434362236536,
            "Po-211": 0.0007424423301461693,
            "Po-214": 84.99649018398776,
            "Po-215": 0.26900084425583065,
            "Po-218": 85.01434319248591,
            "Ra-223": 0.26900062820528614,
            "Ra-226": 85.01434319228659,
            "Rn-218": 1.7002868638497185e-05,
            "Rn-219": 0.26900062820528614,
            "Rn-222": 85.0143431924858,
            "Th-227": 0.2652884195245263,
            "Th-230": 85.01431274847525,
            "Th-231": 0.26898810215560653,
            "Th-234": 85.00820732310407,
            "Tl-206": 0.00011383420610068998,
            "Tl-207": 0.26825840192571576,
            "Tl-210": 0.01785300849981999,
            "U-234": 85.01287846492669,
            "U-235": 0.2689881021544942,
            "U-238": 85.00820732184867,
        }
        if calculated != expected:
            warnings.warn(warning_message_if_dict_not_equal(calculated, expected))
        dict_assert_almost_equal(self, calculated, expected)

    def test_cumulative_decays(self) -> None:
        """
        Test Inventory.cumulative_decays() calculations.
        """

        inv = Inventory({"H-3": 10.0}, "num")
        self.assertEqual(inv.cumulative_decays(12.32, "y"), {"H-3": 5.0})
        self.assertEqual(inv.cumulative_decays(1e6, "y"), {"H-3": 10.0})

        inv = Inventory({"Sr-90": 10.0}, "num")
        self.assertEqual(
            inv.cumulative_decays(1e6, "y"), {"Sr-90": 10.0, "Y-90": 10.000000000000002}
        )

        inv = Inventory({"Tc-99m": 2.3, "I-123": 5.8}, "num")
        calculated = inv.cumulative_decays(20.0, "h")
        expected = {
            "I-123": 3.759540755465226,
            "Te-123m": 4.72801274418656e-07,
            "Te-123": -9.485843969524933e-18,
            "Tc-99m": 2.0704925198993647,
            "Tc-99": 1.0500063924036559e-08,
        }
        if calculated != expected:
            warnings.warn(warning_message_if_dict_not_equal(calculated, expected))
        dict_assert_almost_equal(self, calculated, expected)

        inv = Inventory({"U-238": 99.274, "U-235": 0.720, "U-234": 0.005}, "num")
        calculated = inv.cumulative_decays(1e9, "y")
        expected = {
            "Ac-227": 0.45099937182594435,
            "Th-227": 0.44477558047547366,
            "Ra-226": 14.264656807713404,
            "Fr-223": 0.006223791331197811,
            "Ra-223": 0.4509993717947137,
            "Rn-222": 14.264656807514221,
            "Rn-219": 0.4509993717947137,
            "At-219": 3.734274798718682e-07,
            "Po-218": 14.264656807514106,
            "At-218": 0.002852931361502821,
            "Rn-218": 2.852931361502822e-06,
            "Bi-215": 3.622246554757077e-07,
            "Po-215": 0.4509997340193692,
            "Pb-214": 14.261803876151637,
            "Bi-214": 14.264653954581055,
            "Po-214": 14.261661230181954,
            "Pb-211": 0.450999734019343,
            "Bi-211": 0.4509997340193414,
            "Po-211": 0.0012447592658933822,
            "Tl-210": 0.0029955773304620116,
            "Pb-210": 14.264656385104534,
            "Bi-210": 14.264656384843379,
            "Po-210": 14.26465637763465,
            "U-238": 14.265792678151335,
            "Tl-207": 0.4497549747534446,
            "Hg-206": 2.710284713169805e-07,
            "Tl-206": 1.910037489931004e-05,
            "U-235": 0.4510118978455057,
            "Th-234": 14.26579267689593,
            "Pa-234m": 14.265792676895886,
            "Pa-234": 0.022825268283010153,
            "U-234": 14.266121535073315,
            "Th-231": 0.45101189784439344,
            "Pa-231": 0.45099938014509444,
            "Th-230": 14.264687251524753,
        }
        if calculated != expected:
            warnings.warn(warning_message_if_dict_not_equal(calculated, expected))
        dict_assert_almost_equal(self, calculated, expected)

    def test_half_lives(self) -> None:
        """
        Test method to fetch half-lives of nuclides in the Inventory.
        """

        inv = Inventory({"C-14": 1.0, "H-3": 2.0}, "num")
        self.assertEqual(inv.half_lives("y"), {"C-14": 5700.0, "H-3": 12.32})
        self.assertEqual(
            inv.half_lives("readable"), {"C-14": "5.70 ky", "H-3": "12.32 y"}
        )

    def test_progeny(self) -> None:
        """
        Test method to fetch progeny of nuclides in the Inventory.
        """

        inv = Inventory({"C-14": 1.0, "K-40": 2.0}, "num")
        self.assertEqual(inv.progeny(), {"C-14": ["N-14"], "K-40": ["Ca-40", "Ar-40"]})

    def test_branching_fractions(self) -> None:
        """
        Test method to fetch branching fractions of nuclides in the Inventory.
        """

        inv = Inventory({"C-14": 1.0, "K-40": 2.0}, "num")
        self.assertEqual(
            inv.branching_fractions(), {"C-14": [1.0], "K-40": [0.8914, 0.1086]}
        )

    def test_decay_modes(self) -> None:
        """
        Test method to fetch decay modes of nuclides in the Inventory.
        """

        inv = Inventory({"C-14": 1.0, "K-40": 2.0}, "num")
        self.assertEqual(
            inv.decay_modes(),
            {"C-14": ["\u03b2-"], "K-40": ["\u03b2-", "\u03b2+ \u0026 EC"]},
        )

    def test_decay_time_series_pandas(self) -> None:
        """
        Test method to fetch data of nuclides in the Inventory
        """
        inv = Inventory({"C-14": 1.0})
        self.assertIsNone(
            pd.testing.assert_frame_equal(
                inv.decay_time_series_pandas(
                    time_period=10, time_units="ky", decay_units="mass_frac", npoints=4
                ),
                pd.DataFrame(
                    data={
                        "C-14": [
                            1.0,
                            0.6667465897861368,
                            0.4445504227269143,
                            0.2964018201597633,
                        ],
                        "N-14": [
                            0.0,
                            0.3332534102138631,
                            0.5554495772730857,
                            0.7035981798402366,
                        ],
                        "Time (ky)": [0.0, 3.3333333333333335, 6.666666666666667, 10.0],
                    }
                ).set_index("Time (ky)"),
            )
        )

    def test_decay_time_series(self) -> None:
        """
        Test method to fetch decay data of nuclides in the Inventory as list and dict tuple
        """
        inv = Inventory({"C-14": 1.0})
        time, data = inv.decay_time_series(
            time_period=10, time_units="ky", decay_units="mass_frac", npoints=4
        )
        self.assertEqual(time, [0.0, 3.3333333333333335, 6.666666666666667, 10.0])
        self.assertEqual(
            data,
            {
                "C-14": [
                    1.0,
                    0.6667465897861368,
                    0.4445504227269143,
                    0.2964018201597633,
                ],
                "N-14": [
                    0.0,
                    0.3332534102138631,
                    0.5554495772730857,
                    0.7035981798402366,
                ],
            },
        )

    @patch("matplotlib.pyplot.show")
    def test_plot(self, mock_show) -> None:
        """
        Test method to create decay plots.
        """

        inv = Inventory({"C-14": 1.0, "K-40": 2.0}, "MBq")
        _, axes = inv.plot(105, "ky", yunits="kBq")
        self.assertEqual(axes.get_xscale(), "linear")
        self.assertEqual(axes.get_yscale(), "linear")
        self.assertEqual(axes.get_xlabel(), "Time (ky)")
        self.assertEqual(axes.get_ylabel(), "Activity (kBq)")
        self.assertEqual(axes.get_xlim(), (-5.25, 110.25))
        self.assertEqual(axes.get_ylim(), (0.0, 2100.0))
        self.assertEqual(
            axes.get_legend_handles_labels()[-1],
            ["K-40", "Ca-40", "Ar-40", "C-14", "N-14"],
        )

        inv = Inventory({"C-14": 1.0, "K-40": 2.0}, "kg")
        _, axes = inv.plot(
            100, "ky", xmin=50, ymin=1.0, ymax=2.5, yunits="kg", display="K40"
        )
        self.assertEqual(axes.get_xlim(), (47.5, 102.5))
        self.assertEqual(axes.get_ylim(), (1.0, 2.5))
        self.assertEqual(axes.get_ylabel(), "Mass (kg)")
        self.assertEqual(axes.get_legend_handles_labels()[-1], ["K-40"])

        inv = Inventory({"C-14": 1.0, "K-40": 2.0}, "num")
        _, axes = inv.plot(100, "ky", yunits="num", order="alphabetical")
        self.assertEqual(
            axes.get_legend_handles_labels()[-1],
            ["Ar-40", "C-14", "Ca-40", "K-40", "N-14"],
        )
        self.assertEqual(axes.get_ylabel(), "Number of atoms")

        inv = Inventory({"C-14": 1.0, "K-40": 2.0}, "Bq")
        _, axes = inv.plot(100, "ky", yunits="activity_frac")
        self.assertEqual(axes.get_ylabel(), "Activity fraction")

        inv = Inventory({"C-14": 1.0, "K-40": 2.0}, "g")
        _, axes = inv.plot(100, "ky", yunits="mass_frac")
        self.assertEqual(axes.get_ylabel(), "Mass fraction")

        inv = Inventory({"C-14": 1.0, "K-40": 2.0}, "kmol")
        _, axes = inv.plot(100, "ky", yunits="mol_frac")
        self.assertEqual(axes.get_ylabel(), "Mole fraction")

        with self.assertRaises(ValueError):
            inv.plot(100, "ky", order="invalid")
        with self.assertRaises(ValueError):
            inv.plot(100, "ky", yunits="invalid")

    @patch("radioactivedecay.inventory._write_csv_file")
    def test_to_csv_invalid_units(self, mock__write_csv_file) -> None:
        """
        Test writing Inventory to csv file with invalid units.
        """

        inv = Inventory({"H-3": 10.0, "14C": 50.0})
        with self.assertRaises(ValueError):
            inv.to_csv("test_file.csv", units="fake")

    @patch("radioactivedecay.inventory._write_csv_file")
    def test_to_csv_write_units(self, mock__write_csv_file) -> None:
        """
        Test writing Inventory to csv file including writing units.
        """

        inv = Inventory({"H-3": 10.0, "14C": 50.0})
        kwargs = {
            "filename": "test_file.csv",
            "delimiter": "\t",
            "write_units": True,
            "encoding": "utf-16",
        }
        inv.to_csv(**kwargs)
        mock__write_csv_file.assert_called_once_with(
            kwargs["filename"],
            [["C-14", "50.0", "Bq"], ["H-3", "10.0", "Bq"]],
            kwargs["delimiter"],
            kwargs["encoding"],
        )

    @patch("radioactivedecay.inventory._write_csv_file")
    def test_to_csv_dont_write_units(self, mock__write_csv_file) -> None:
        """
        Test writing Inventory to csv file without writing units.
        """

        inv = Inventory({"H-3": 10.0}, "kg")
        kwargs = {
            "filename": "test_file.csv",
            "units": "g",
        }
        inv.to_csv(**kwargs)
        mock__write_csv_file.assert_called_with(
            kwargs["filename"],
            [["H-3", "10000.0"]],
            ",",
            "utf-8",
        )

        inv = Inventory({"H-3": 20.2}, "mol")
        kwargs = {
            "filename": "test_file.csv",
            "units": "mmol",
        }
        inv.to_csv(**kwargs)
        mock__write_csv_file.assert_called_with(
            kwargs["filename"],
            [["H-3", "20200.0"]],
            ",",
            "utf-8",
        )

        inv = Inventory({"H-3": 20.2}, "num")
        kwargs = {
            "filename": "test_file.csv",
            "units": "num",
        }
        inv.to_csv(**kwargs)
        mock__write_csv_file.assert_called_with(
            kwargs["filename"],
            [["H-3", "20.2"]],
            ",",
            "utf-8",
        )

    def test___eq__(self) -> None:
        """
        Test Inventory equality.
        """

        inv1 = Inventory({"H-3": 10.0})
        inv2 = Inventory({"H3": 10.0})
        self.assertEqual(inv1, inv2)

        decay_data = load_dataset("icrp107_ame2020_nubase2020", load_sympy=True)
        inv2 = Inventory({"H-3": 10.0}, decay_data=decay_data)
        self.assertEqual(inv1, inv2)

        self.assertFalse(inv1 == "random object")

    def test___ne__(self) -> None:
        """
        Test Inventory inequality.
        """

        inv1 = Inventory({"H-3": 10.0})
        inv2 = Inventory({"Cs-137": 10.0})
        self.assertNotEqual(inv1, inv2)

        inv1 = Inventory({"H-3": 10.0})
        inv2 = Inventory({"H-3": 5.0})
        self.assertNotEqual(inv1, inv2)

        self.assertTrue(inv1 != "random object")

    def test___repr__(self) -> None:
        """
        Test Inventory __repr__ strings.
        """

        inv = Inventory({"H-3": 10.0}, "Bq")
        self.assertEqual(
            repr(inv),
            "Inventory activities (Bq): {'H-3': 10.0}, decay dataset: icrp107_ame2020_nubase2020",
        )


class TestInventoryHP(unittest.TestCase):
    """
    Unit tests for the inventory.py InventoryHP class.
    """

    def test_instantiation(self) -> None:
        """
        Test InventoryHP instantiation.
        """

        inv = InventoryHP({"H3": 1}, "Bq", True)
        self.assertEqual(
            inv.contents, {"H-3": Integer(242988330816) / (Integer(625) * log(2))}
        )
        self.assertEqual(inv.decay_data, DEFAULTDATA)
        self.assertEqual(inv.decay_matrices, DEFAULTDATA.sympy_data)

        temp_data = copy.deepcopy(DEFAULTDATA)
        backup_year_conv = temp_data._sympy_year_conv
        temp_data._sympy_year_conv = None
        with self.assertRaises(ValueError):
            InventoryHP({"H3": 1}, "Bq", True, temp_data)
        temp_data._sympy_year_conv = backup_year_conv
        temp_data._sympy_data = None
        with self.assertRaises(ValueError):
            InventoryHP({"H3": 1}, "Bq", True, temp_data)

        # check instantiation with a stable nuclide activity raises a ValueError
        with self.assertRaises(ValueError):
            InventoryHP({"He-3": 0.0}, "Bq", False)

    def test_numbers(self) -> None:
        """
        Test InventoryHP.numbers() method.
        """

        inv = InventoryHP({"H-3": 1}, "num")
        self.assertEqual(inv.numbers(), {"H-3": 1.0})
        inv = InventoryHP({"Tc-99m": 2.3, "I-123": 5.8})
        self.assertEqual(
            inv.numbers(), {"I-123": 399738.47946141585, "Tc-99m": 71852.27235544211}
        )

    def test_activities(self) -> None:
        """
        Test InventoryHP.activities() method.
        """

        inv = InventoryHP({"H-3": 1})
        self.assertEqual(inv.activities(), {"H-3": 1})
        inv = InventoryHP({"Tc-99m": 2.3, "I-123": 5.8})
        self.assertEqual(inv.activities(), {"I-123": 5.8, "Tc-99m": 2.3})
        inv = InventoryHP({"Tc-99m": 2300, "I-123": 5800})
        self.assertEqual(inv.activities("kBq"), {"I-123": 5.8, "Tc-99m": 2.3})

    def test_masses(self) -> None:
        """
        Test InventoryHP.masses() method.
        """

        inv = InventoryHP({"H-3": 1}, "g")
        self.assertEqual(inv.masses(), {"H-3": 1})
        inv = InventoryHP({"Tc-99m": 2.3, "I-123": 5.8})
        self.assertEqual(
            inv.masses(),
            {"I-123": 8.158243973887584e-17, "Tc-99m": 1.1800869622748503e-17},
        )
        self.assertEqual(
            inv.masses("pg"),
            {"I-123": 8.158243973887584e-05, "Tc-99m": 1.1800869622748502e-05},
        )

    def test_decay(self) -> None:
        """
        Test InventoryHP.decay() calculations.
        """

        inv = InventoryHP({"U-238": 99.274, "U-235": 0.720, "U-234": 0.005}, "Bq")
        self.assertEqual(
            inv.decay(1e9, "y").activities(),
            {
                "Ac-227": 0.26900062817405557,
                "At-218": 0.01700286863849718,
                "At-219": 2.227325201281318e-07,
                "Bi-210": 85.01434361515662,
                "Bi-211": 0.2690008442558584,
                "Bi-214": 85.01432618961894,
                "Bi-215": 2.1605054452429227e-07,
                "Fr-223": 0.003712208668802187,
                "Hg-206": 1.6152725286830195e-06,
                "Pa-231": 0.2690006198549054,
                "Pa-234": 0.13601313171698984,
                "Pa-234m": 85.00820732310412,
                "Pb-206": 0.0,
                "Pb-207": 0.0,
                "Pb-210": 85.01434361489547,
                "Pb-211": 0.26900084425585685,
                "Pb-214": 84.99734032384836,
                "Po-210": 85.01434362236536,
                "Po-211": 0.0007424423301461693,
                "Po-214": 84.99649018398776,
                "Po-215": 0.26900084425583065,
                "Po-218": 85.0143431924859,
                "Ra-223": 0.2690006282052861,
                "Ra-226": 85.0143431922866,
                "Rn-218": 1.7002868638497178e-05,
                "Rn-219": 0.26900062820528614,
                "Rn-222": 85.01434319248578,
                "Th-227": 0.26528841952452625,
                "Th-230": 85.01431274847525,
                "Th-231": 0.26898810215560653,
                "Th-234": 85.00820732310407,
                "Tl-206": 0.00011383420610068996,
                "Tl-207": 0.2682584019257157,
                "Tl-210": 0.017853008499819988,
                "U-234": 85.01287846492669,
                "U-235": 0.26898810215449415,
                "U-238": 85.00820732184867,
            },
        )

        inv.sig_fig = 0
        with self.assertRaises(ValueError):
            inv.decay(1e9, "y")

    def test_cumulative_decays(self) -> None:
        """
        Test InventoryHP.cumulative_decays() calculations.
        """

        inv = InventoryHP({"H-3": 10.0}, "num")
        self.assertEqual(inv.cumulative_decays(12.32, "y"), {"H-3": 5.0})
        self.assertEqual(inv.cumulative_decays(1e6, "y"), {"H-3": 10.0})

        inv = InventoryHP({"Sr-90": 10.0}, "num")
        self.assertEqual(inv.cumulative_decays(1e6, "y"), {"Sr-90": 10.0, "Y-90": 10.0})

        inv = InventoryHP({"Tc-99m": 2.3, "I-123": 5.8}, "num")
        self.assertEqual(
            inv.cumulative_decays(20.0, "h"),
            {
                "I-123": 3.7595407554652263,
                "Te-123m": 4.728012744186443e-07,
                "Te-123": 5.801848145646267e-18,
                "Tc-99m": 2.0704925198993647,
                "Tc-99": 1.0500063848957902e-08,
            },
        )

        inv = InventoryHP({"U-238": 99.274, "U-235": 0.720, "U-234": 0.005}, "num")
        self.assertEqual(
            inv.cumulative_decays(1e9, "y"),
            {
                "Ac-227": 0.45099937182594446,
                "Th-227": 0.4447755804754738,
                "Ra-226": 14.264656807713404,
                "Fr-223": 0.0062237913311978125,
                "Ra-223": 0.45099937179471394,
                "Rn-222": 14.264656807514218,
                "Rn-219": 0.4509993717947139,
                "At-219": 3.734274798718682e-07,
                "Po-218": 14.264656807514106,
                "At-218": 0.002852931361502821,
                "Rn-218": 2.8529313615028208e-06,
                "Bi-215": 3.6222465547570775e-07,
                "Po-215": 0.4509997340193693,
                "Pb-214": 14.261803876151633,
                "Bi-214": 14.264653954581055,
                "Po-214": 14.261661230181954,
                "Pb-211": 0.4509997340193431,
                "Bi-211": 0.45099973401934157,
                "Po-211": 0.0012447592658933826,
                "Tl-210": 0.0029955773304620116,
                "Pb-210": 14.264656385104532,
                "Bi-210": 14.264656384843379,
                "Po-210": 14.264656377634646,
                "U-238": 14.265792678151334,
                "Tl-207": 0.4497549747534447,
                "Hg-206": 2.710284713169805e-07,
                "Tl-206": 1.9100374899310035e-05,
                "U-235": 0.4510118978455059,
                "Th-234": 14.265792676895929,
                "Pa-234m": 14.265792676895886,
                "Pa-234": 0.02282526828301015,
                "U-234": 14.266121535073315,
                "Th-231": 0.4510118978443935,
                "Pa-231": 0.45099938014509455,
                "Th-230": 14.264687251524753,
            },
        )

        # Catch incorrect sig_fig or no SymPy data in decay dataset
        inv.sig_fig = 0
        with self.assertRaises(ValueError):
            inv.cumulative_decays(1e9, "y")

    @patch("matplotlib.pyplot.show")
    def test_plot(self, mock_show) -> None:
        """
        Test InventoryHP.plot() method.
        """

        inv = InventoryHP({"C-14": 1.0, "K-40": 2.0}, "mol")
        _, axes = inv.plot(
            100,
            xscale="log",
            yscale="log",
            yunits="mmol",
            display=["K40", "C14"],
        )
        self.assertEqual(axes.get_xscale(), "log")
        self.assertEqual(axes.get_yscale(), "log")
        self.assertEqual(axes.get_xlabel(), "Time (s)")
        self.assertEqual(axes.get_ylabel(), "Number of moles (mmol)")
        self.assertEqual(axes.get_xlim()[0], 0.0707945784384138)
        self.assertEqual(axes.get_ylim(), (949.999999633917, 2100.0))
        self.assertEqual(axes.get_legend_handles_labels()[-1], ["K-40", "C-14"])

    @patch("radioactivedecay.inventory._write_csv_file")
    def test_to_csv_write_units(self, mock__write_csv_file) -> None:
        """
        Test writing Inventory to csv file including writing units.
        """

        inv = InventoryHP({"H-3": 10.0, "14C": 50.0})
        kwargs = {
            "filename": "test_file.csv",
            "delimiter": "\t",
            "write_units": True,
            "header": ["nuclide", "quantity", "units"],
            "encoding": "utf-16",
        }
        inv.to_csv(**kwargs)
        mock__write_csv_file.assert_called_once_with(
            kwargs["filename"],
            [
                ["nuclide", "quantity", "units"],
                ["C-14", "50.0", "Bq"],
                ["H-3", "10.0", "Bq"],
            ],
            kwargs["delimiter"],
            kwargs["encoding"],
        )

    def test___repr__(self) -> None:
        """
        Test InventoryHP __repr__ strings.
        """

        inv = InventoryHP({"H-3": 10.0}, "Bq")
        self.assertEqual(
            repr(inv),
            "InventoryHP activities (Bq): {'H-3': 10.0}, decay dataset: icrp107_ame2020_nubase2020",
        )


if __name__ == "__main__":
    unittest.main()
