"""Utility functions to be used in plotting."""

import re
import numpy as np


def axis_label_in_latex(label_text, unit, only_text=True):
    """
    Get axis label for plotly plots that can show units in latex.

    Parameters
    ----------
    label_text : str
        Text to show on label, may be expressed in latex
    unit : astropy.units
        Unit of the label which needs to be expressed in latex
    only_text : bool
        If label_text is expressed purely in text (i.e. without
        using latex) or not. Default value is True

    Returns
    -------
    str
        Latex string for label renderable by plotly
    """
    unit_in_latex = unit.to_string("latex_inline").strip("$")

    # If present, place s^{-1} just after erg
    if "erg" in unit_in_latex and "s^{-1}" in unit_in_latex:
        constituent_units = (
            re.compile("\\\mathrm\{(.*)\}")
            .findall(unit_in_latex)[0]
            .split("\\,")
        )
        constituent_units.remove("s^{-1}")
        constituent_units.insert(constituent_units.index("erg") + 1, "s^{-1}")
        constituent_units_string = "\\,".join(constituent_units)
        unit_in_latex = f"\\mathrm{{{constituent_units_string}}}"

    if only_text:
        return f"$\\text{{{label_text}}}\\,[{unit_in_latex}]$"
    else:
        return f"${label_text}\\,[{unit_in_latex}]$"


def get_mid_point_idx(arr):
    """
    Get index of the middle point of a sorted array (ascending or descending).

    The values in array may not be evenly distributed so it picks the middle
    point not by index but by their values.

    Parameters
    ----------
    arr : np.array

    Returns
    -------
    int
    """
    mid_value = (arr[0] + arr[-1]) / 2
    return np.abs(arr - mid_value).argmin()


def to_rgb255_string(color_tuple):
    """
    Convert a matplotlib RGBA tuple to a generic RGB 255 string.

    Parameters
    ----------
    color_tuple : tuple
        Matplotlib RGBA tuple of float values in closed interval [0, 1]

    Returns
    -------
    str
        RGB string of format rgb(r,g,b) where r,g,b are integers between
        0 and 255 (both inclusive)
    """
    color_tuple_255 = tuple([int(x * 255) for x in color_tuple[:3]])
    return f"rgb{color_tuple_255}"
