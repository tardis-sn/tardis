"""
Unit tests for nuclide.py functions, classes and methods.
"""

import unittest

from radioactivedecay.decaydata import load_dataset
from radioactivedecay.nuclide import Nuclide


class TestNuclide(unittest.TestCase):
    """
    Unit tests for the nuclide.py Nuclide class.
    """

    def test_nuclide_instantiation(self) -> None:
        """
        Test instantiation of Nuclide objects.
        """
        decay_data = load_dataset("icrp107_ame2020_nubase2020", load_sympy=True)

        nuc = Nuclide("Rn-222")
        self.assertEqual(nuc.nuclide, "Rn-222")
        self.assertEqual(nuc.decay_data, decay_data)

        nuc = Nuclide("222Rn")
        self.assertEqual(nuc.nuclide, "Rn-222")

        nuc = Nuclide(611450000)
        self.assertEqual(nuc.nuclide, "Pm-145")

    def test_nuclide_Z(self) -> None:
        """
        Test Nuclide Z property.
        """

        nuc = Nuclide("H-3")
        self.assertEqual(nuc.Z, 1)

    def test_nuclide_A(self) -> None:
        """
        Test Nuclide A property.
        """

        nuc = Nuclide("H-3")
        self.assertEqual(nuc.A, 3)

    def test_nuclide_state(self) -> None:
        """
        Test Nuclide state property.
        """

        nuc = Nuclide("H-3")
        self.assertEqual(nuc.state, "")

        nuc = Nuclide("Ba-137m")
        self.assertEqual(nuc.state, "m")

        nuc = Nuclide("Ir-192n")
        self.assertEqual(nuc.state, "n")

    def test_nuclide_id(self) -> None:
        """
        Test Nuclide id property.
        """

        nuc = Nuclide("H-3")
        self.assertEqual(nuc.id, 10030000)

        nuc = Nuclide("I-118m")
        self.assertEqual(nuc.id, 531180001)

        nuc = Nuclide(190400000)
        self.assertEqual(nuc.id, 190400000)

    def test_nuclide_atomic_mass(self) -> None:
        """
        Test Nuclide atomic_mass property.
        """

        nuc = Nuclide("H-3")
        self.assertEqual(nuc.atomic_mass, 3.01604928132)

    def test_nuclide_half_life(self) -> None:
        """
        Test Nuclide half_life() method.
        """

        nuc = Nuclide("H-3")
        self.assertEqual(nuc.half_life(), 388781329.30560005)
        self.assertEqual(nuc.half_life("y"), 12.32)
        self.assertEqual(nuc.half_life("readable"), "12.32 y")

    def test_nuclide_progeny(self) -> None:
        """
        Test Nuclide progeny() method.
        """

        nuc = Nuclide("K-40")
        self.assertEqual(nuc.progeny()[0], "Ca-40")
        self.assertEqual(nuc.progeny()[1], "Ar-40")

    def test_nuclide_branching_fractions(self) -> None:
        """
        Test Nuclide branching_fractions() method.
        """

        nuc = Nuclide("K-40")
        self.assertEqual(nuc.branching_fractions()[0], 0.8914)
        self.assertEqual(nuc.branching_fractions()[1], 0.1086)

    def test_nuclide_decay_modes(self) -> None:
        """
        Test Nuclide decay_modes() method.
        """

        nuc = Nuclide("K-40")
        self.assertEqual(nuc.decay_modes()[0], "\u03b2-")
        self.assertEqual(nuc.decay_modes()[1], "\u03b2+ & EC")

    def test_nuclide_plot(self) -> None:
        """
        Test Nuclide plot() method.

        Only testing auto-generation of limits so far.
        """

        nuc = Nuclide("H-3")
        _, axes = nuc.plot()
        self.assertEqual(axes.get_xlim(), (-0.3, 0.3))
        self.assertEqual(axes.get_ylim(), (-1.3, 0.3))

        nuc = Nuclide("Mo-99")
        _, axes = nuc.plot()
        self.assertEqual(axes.get_xlim(), (-0.3, 1.3))
        self.assertEqual(axes.get_ylim(), (-2.3, 0.3))

        nuc = Nuclide("Es-256")
        _, axes = nuc.plot()
        self.assertEqual(axes.get_xlim(), (-0.3, 2.3))
        self.assertEqual(axes.get_ylim(), (-19.3, 0.3))

        nuc = Nuclide("Cu-64")
        _, axes = nuc.plot()
        self.assertEqual(axes.get_xlim(), (-0.3, 1.3))
        self.assertEqual(axes.get_ylim(), (-1.3, 0.3))

    def test_nuclide___repr__(self) -> None:
        """
        Test Nuclide __repr__ strings.
        """

        nuc = Nuclide("H-3")
        self.assertEqual(
            repr(nuc),
            "Nuclide: H-3, decay dataset: icrp107_ame2020_nubase2020",
        )

    def test_nuclide___eq__(self) -> None:
        """
        Test Nuclide equality.
        """

        nuc1 = Nuclide("K-40")
        nuc2 = Nuclide("40K")
        nuc3 = Nuclide(190400000)
        self.assertEqual(nuc1, nuc2)
        self.assertEqual(nuc1, nuc3)

        decay_data = load_dataset("icrp107_ame2020_nubase2020", load_sympy=True)
        nuc2 = Nuclide("K-40", decay_data)
        self.assertEqual(nuc1, nuc2)

        self.assertFalse(nuc1 == "random object")

    def test_nuclide___ne__(self) -> None:
        """
        Test Nuclide inequality.
        """

        nuc1 = Nuclide("K-40")
        nuc2 = Nuclide("H-3")
        self.assertNotEqual(nuc1, nuc2)

        self.assertTrue(nuc1 != "random object")

    def test_nuclide___hash__(self) -> None:
        """
        Test Nuclide hash function.
        """

        nuc = Nuclide("K-40")
        decay_data = load_dataset("icrp107_ame2020_nubase2020", load_sympy=True)
        self.assertEqual(hash(nuc), hash(("K-40", decay_data.dataset_name)))


if __name__ == "__main__":
    unittest.main()
