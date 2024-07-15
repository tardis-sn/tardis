"""
Unit tests for fileio.py functions.
"""

import unittest
from unittest.mock import patch

from radioactivedecay.decaydata import load_dataset
from radioactivedecay.fileio import read_csv
from radioactivedecay.inventory import Inventory, InventoryHP


class TestReadCSV(unittest.TestCase):
    """
    Unit tests for the read_csv() function.
    """

    @patch("radioactivedecay.fileio._read_csv_file")
    def test_inventory(self, mock__read_csv_file) -> None:
        """Test normal inventory from file."""

        mock__read_csv_file.return_value = [["H-3", "1"]]

        inv = read_csv("fake_file.csv")
        self.assertEqual(inv, Inventory({"H-3": 1.0}))

    @patch("radioactivedecay.fileio._read_csv_file")
    def test_multiple_rows(self, mock__read_csv_file) -> None:
        """Test file with more than one nuclide."""

        mock__read_csv_file.return_value = [["H-3", "1"], ["C-14", "2"]]

        inv = read_csv("fake_file.csv")
        self.assertEqual(inv, Inventory({"H-3": 1.0, "C-14": 2.0}))

    @patch("radioactivedecay.fileio._read_csv_file")
    def test_inventory_canonical_id(self, mock__read_csv_file) -> None:
        """Test file with canonical id of nuclide."""
        mock__read_csv_file.return_value = [["10030000", "1"]]

        inv = read_csv("fake_file.csv")
        self.assertEqual(inv, Inventory({"H-3": 1.0}))

    @patch("radioactivedecay.fileio._read_csv_file")
    def test_invalid_rows(self, mock__read_csv_file) -> None:
        "Test file with invalid rows."

        mock__read_csv_file.return_value = [["H-31"]]
        with self.assertRaises(ValueError):
            read_csv("fake_file.csv")

        mock__read_csv_file.return_value = [["H-3,1,Bq,mol"]]
        with self.assertRaises(ValueError):
            read_csv("fake_file.csv")

    @patch("radioactivedecay.fileio._read_csv_file")
    def test_empty_file(self, mock__read_csv_file) -> None:
        """Test catches empty file input."""

        mock__read_csv_file.return_value = []
        with self.assertRaises(ValueError):
            read_csv("fake_file.csv")

    @patch("radioactivedecay.fileio._read_csv_file")
    def test_inventoryhp(self, mock__read_csv_file) -> None:
        """Test high precision inventory from file."""

        mock__read_csv_file.return_value = [["H-3", "1"]]

        inv = read_csv("fake_file.csv", "InventoryHP")
        self.assertEqual(inv, InventoryHP({"H-3": 1.0}))

    @patch("radioactivedecay.fileio._read_csv_file")
    def test_invalid_inventory_parameter(self, mock__read_csv_file) -> None:
        """Test invalid inventory parameter."""

        mock__read_csv_file.return_value = [["H-3", "1"]]

        with self.assertRaises(ValueError):
            read_csv("fake_file.csv", "random string")

    @patch("radioactivedecay.fileio._read_csv_file")
    def test_units_parameter(self, mock__read_csv_file) -> None:
        """Test passing units as parameter."""

        mock__read_csv_file.return_value = [["H-3", "1"]]

        inv = read_csv("fake_file.csv", units="mol")
        self.assertEqual(inv, Inventory({"H-3": 1.0}, "mol"))

    @patch("radioactivedecay.fileio._read_csv_file")
    def test_units_in_file(self, mock__read_csv_file) -> None:
        """Test reading units from row in file."""

        mock__read_csv_file.return_value = [["H-3", "1", "mol"]]

        inv = read_csv("fake_file.csv", units="mol")
        self.assertEqual(inv, Inventory({"H-3": 1.0}, "mol"))

    @patch("radioactivedecay.fileio._read_csv_file")
    def test_different_units_in_file(self, mock__read_csv_file) -> None:
        """Test file with multiple nuclides, different units."""

        mock__read_csv_file.return_value = [["H-3", "1", "Ci"], ["C-14", "2", "Bq"]]

        inv = read_csv("fake_file.csv", units="mol")
        self.assertEqual(inv, Inventory({"H-3": 3.7e10, "C-14": 2.0}, "Bq"))

    @patch("radioactivedecay.fileio._read_csv_file")
    def test_units_in_file_overrides(self, mock__read_csv_file) -> None:
        """Test units in file overrides parameter set value."""

        mock__read_csv_file.return_value = [["H-3", "1", "mol"]]

        inv = read_csv("fake_file.csv", units="Ci")
        self.assertEqual(inv, Inventory({"H-3": 1.0}, "mol"))

    @patch("radioactivedecay.fileio._read_csv_file")
    def test_decay_data_parameter(self, mock__read_csv_file) -> None:
        """Test passing DecayData instance as parameter."""

        mock__read_csv_file.return_value = [["H-3", "1"]]

        decay_data = load_dataset("icrp107_ame2020_nubase2020", load_sympy=True)
        inv = read_csv("fake_file.csv", decay_data=decay_data)
        self.assertEqual(inv, Inventory({"H-3": 1.0}, decay_data=decay_data))

    @patch("radioactivedecay.fileio._read_csv_file")
    def test_skiprows_parameter(self, mock__read_csv_file) -> None:
        """Test skipping rows at head of file."""

        mock__read_csv_file.return_value = [["header"], ["H-3", "1"]]

        inv = read_csv("fake_file.csv", skip_rows=1)
        self.assertEqual(inv, Inventory({"H-3": 1.0}))
