"""
Unit tests for plots.py functions.
"""

import unittest

import matplotlib
import matplotlib.pyplot as plt

from radioactivedecay.plots import (
    _check_fig_axes,
    _parse_decay_mode_label,
    _parse_nuclide_label,
)


class TestFunctions(unittest.TestCase):
    """
    Unit tests for plots.py functions.
    """

    def test__parse_nuclide_label(self) -> None:
        """
        Test the parsing of nuclide strings for node labels.
        """

        self.assertEqual(_parse_nuclide_label("H-3"), "³H")
        self.assertEqual(_parse_nuclide_label("Be-7"), "⁷Be")
        self.assertEqual(_parse_nuclide_label("C-10"), "¹⁰C")
        self.assertEqual(_parse_nuclide_label("Ne-19"), "¹⁹Ne")
        self.assertEqual(_parse_nuclide_label("I-118"), "¹¹⁸I")
        self.assertEqual(_parse_nuclide_label("Pd-100"), "¹⁰⁰Pd")
        self.assertEqual(_parse_nuclide_label("Cl-34m"), "³⁴ᵐCl")
        self.assertEqual(_parse_nuclide_label("I-118m"), "¹¹⁸ᵐI")
        self.assertEqual(_parse_nuclide_label("Tb-156m"), "¹⁵⁶ᵐTb")
        self.assertEqual(_parse_nuclide_label("Tb-156n"), "¹⁵⁶ⁿTb")
        self.assertEqual(_parse_nuclide_label("In-129p"), "¹²⁹ᵖIn")
        self.assertEqual(_parse_nuclide_label("Lu-177q"), "¹⁷⁷qLu")
        self.assertEqual(_parse_nuclide_label("Lu-177r"), "¹⁷⁷ʳLu")
        self.assertEqual(_parse_nuclide_label("Lu-174x"), "¹⁷⁴ˣLu")
        self.assertEqual(_parse_nuclide_label("SF"), "various")

    def test__parse_decay_mode_label(self) -> None:
        """
        Test the parsing of decay mode strings for edge labels.
        """

        self.assertEqual(_parse_decay_mode_label("α"), "α")
        self.assertEqual(_parse_decay_mode_label("β+"), "β⁺")
        self.assertEqual(_parse_decay_mode_label("2β+"), "2β⁺")
        self.assertEqual(_parse_decay_mode_label("β+p"), "β⁺p")
        self.assertEqual(_parse_decay_mode_label("β+2p"), "β⁺2p")
        self.assertEqual(_parse_decay_mode_label("β+3p"), "β⁺3p")
        self.assertEqual(_parse_decay_mode_label("β+α"), "β⁺α")
        self.assertEqual(_parse_decay_mode_label("β+SF"), "β⁺SF")
        self.assertEqual(_parse_decay_mode_label("β+ & EC"), "β⁺ & EC")
        self.assertEqual(_parse_decay_mode_label("β-"), "β⁻")
        self.assertEqual(_parse_decay_mode_label("2β-"), "2β⁻")
        self.assertEqual(_parse_decay_mode_label("β-n"), "β⁻n")
        self.assertEqual(_parse_decay_mode_label("β-2n"), "β⁻2n")
        self.assertEqual(_parse_decay_mode_label("β-3n"), "β⁻3n")
        self.assertEqual(_parse_decay_mode_label("β-α"), "β⁻α")
        self.assertEqual(_parse_decay_mode_label("β-d"), "β⁻d")
        self.assertEqual(_parse_decay_mode_label("β-t"), "β⁻t")
        self.assertEqual(_parse_decay_mode_label("β-SF"), "β⁻SF")
        self.assertEqual(_parse_decay_mode_label("EC"), "EC")
        self.assertEqual(_parse_decay_mode_label("ε"), "ε")
        self.assertEqual(_parse_decay_mode_label("e+"), "e⁺")
        self.assertEqual(_parse_decay_mode_label("IT"), "IT")
        self.assertEqual(_parse_decay_mode_label("SF"), "SF")
        self.assertEqual(_parse_decay_mode_label("p"), "p")
        self.assertEqual(_parse_decay_mode_label("2p"), "2p")
        self.assertEqual(_parse_decay_mode_label("n"), "n")
        self.assertEqual(_parse_decay_mode_label("2n"), "2n")
        self.assertEqual(_parse_decay_mode_label("12C"), "¹²C")
        self.assertEqual(_parse_decay_mode_label("14C"), "¹⁴C")
        self.assertEqual(_parse_decay_mode_label("20O"), "²⁰O")
        self.assertEqual(_parse_decay_mode_label("23F"), "²³F")
        self.assertEqual(_parse_decay_mode_label("22Ne"), "²²Ne")
        self.assertEqual(_parse_decay_mode_label("24Ne"), "²⁴Ne")
        self.assertEqual(_parse_decay_mode_label("25Ne"), "²⁵Ne")
        self.assertEqual(_parse_decay_mode_label("26Ne"), "²⁶Ne")
        self.assertEqual(_parse_decay_mode_label("28Mg"), "²⁸Mg")
        self.assertEqual(_parse_decay_mode_label("29Mg"), "²⁹Mg")
        self.assertEqual(_parse_decay_mode_label("30Mg"), "³⁰Mg")
        self.assertEqual(_parse_decay_mode_label("32Si"), "³²Si")
        self.assertEqual(_parse_decay_mode_label("34Si"), "³⁴Si")

    def test__check_fig_axes(self) -> None:
        """
        Test the parsing of user-defined Matplotlib Figure and Axes objects.
        """

        fig_in, axes_in = plt.subplots()
        fig, axes = _check_fig_axes(fig_in, axes_in)
        self.assertIsInstance(fig, matplotlib.figure.Figure)
        self.assertIsInstance(axes, matplotlib.axes.Axes)

        fig, axes = _check_fig_axes(fig_in, None)
        self.assertIsInstance(fig, matplotlib.figure.Figure)
        self.assertIsInstance(axes, matplotlib.axes.Axes)

        fig, axes = _check_fig_axes(None, axes_in)
        self.assertIsInstance(fig, matplotlib.figure.Figure)
        self.assertIsInstance(axes, matplotlib.axes.Axes)

        fig, axes = _check_fig_axes(None, None)
        self.assertIsInstance(fig, matplotlib.figure.Figure)
        self.assertIsInstance(axes, matplotlib.axes.Axes)


if __name__ == "__main__":
    unittest.main()
