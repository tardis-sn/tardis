"""
The utils module contains functions to parse nuclide strings and manipulate lists and dictionaries.

The docstring code examples assume that ``radioactivedecay`` has been imported
as ``rd``:

.. highlight:: python
.. code-block:: python

    >>> import radioactivedecay as rd

"""

from typing import Dict, List, Tuple, Union

import numpy as np
from sympy.core.expr import Expr

Z_DICT = {
    1: "H",
    2: "He",
    3: "Li",
    4: "Be",
    5: "B",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
    10: "Ne",
    11: "Na",
    12: "Mg",
    13: "Al",
    14: "Si",
    15: "P",
    16: "S",
    17: "Cl",
    18: "Ar",
    19: "K",
    20: "Ca",
    21: "Sc",
    22: "Ti",
    23: "V",
    24: "Cr",
    25: "Mn",
    26: "Fe",
    27: "Co",
    28: "Ni",
    29: "Cu",
    30: "Zn",
    31: "Ga",
    32: "Ge",
    33: "As",
    34: "Se",
    35: "Br",
    36: "Kr",
    37: "Rb",
    38: "Sr",
    39: "Y",
    40: "Zr",
    41: "Nb",
    42: "Mo",
    43: "Tc",
    44: "Ru",
    45: "Rh",
    46: "Pd",
    47: "Ag",
    48: "Cd",
    49: "In",
    50: "Sn",
    51: "Sb",
    52: "Te",
    53: "I",
    54: "Xe",
    55: "Cs",
    56: "Ba",
    57: "La",
    58: "Ce",
    59: "Pr",
    60: "Nd",
    61: "Pm",
    62: "Sm",
    63: "Eu",
    64: "Gd",
    65: "Tb",
    66: "Dy",
    67: "Ho",
    68: "Er",
    69: "Tm",
    70: "Yb",
    71: "Lu",
    72: "Hf",
    73: "Ta",
    74: "W",
    75: "Re",
    76: "Os",
    77: "Ir",
    78: "Pt",
    79: "Au",
    80: "Hg",
    81: "Tl",
    82: "Pb",
    83: "Bi",
    84: "Po",
    85: "At",
    86: "Rn",
    87: "Fr",
    88: "Ra",
    89: "Ac",
    90: "Th",
    91: "Pa",
    92: "U",
    93: "Np",
    94: "Pu",
    95: "Am",
    96: "Cm",
    97: "Bk",
    98: "Cf",
    99: "Es",
    100: "Fm",
    101: "Md",
    102: "No",
    103: "Lr",
    104: "Rf",
    105: "Db",
    106: "Sg",
    107: "Bh",
    108: "Hs",
    109: "Mt",
    110: "Ds",
    111: "Rg",
    112: "Cn",
    113: "Nh",
    114: "Fl",
    115: "Mc",
    116: "Lv",
    117: "Ts",
    118: "Og",
}
SYM_DICT = dict((v, k) for k, v in Z_DICT.items())
METASTABLE_CHARS = ["m", "n", "p", "q", "r", "x"]


def get_metastable_chars() -> List[str]:
    """
    Returns list of allowed metastable state characters. Currently up to sixth metastable state is
    supported, based on the existence of sixth metastable state isomers in NUBASE2020, e.g.
    Lu-174x.

    Note metastable state chars are "m", "n", "p", "q", "r", "x" for the first to sixth metastable
    state, respectively, i.e. "o" is not a metastable state char, nor letters between "r" and "x".
    See NUBASE2020 paper (F.G. Kondev et al 2021 Chinese Phys. C 45 030001) for more details.
    """

    return METASTABLE_CHARS


def Z_to_elem(Z: int) -> str:
    """
    Converts atomic number to element symbol.

    Parameters
    ----------
    Z : int
        Atomic number.

    Returns
    -------
    str
        Element string.

    Examples
    --------
    >>> rd.utils.Z_to_elem(1)
    'H'
    >>> rd.utils.Z_to_elem(35)
    'Br'

    """

    return Z_DICT[Z]


def elem_to_Z(sym: str) -> int:
    """
    Converts element symbol to atomic number.

    Parameters
    ----------
    sym : str
        Element string.

    Returns
    -------
    int
        Atomic number.

    Examples
    --------
    >>> rd.utils.elem_to_Z('H')
    1
    >>> rd.utils.elem_to_Z('Br')
    35

    """

    return SYM_DICT[sym]


def build_id(Z: int, A: int, state: str = "") -> int:
    """
    Builds a canonical nuclide id from atomic number, atomic mass, and
    and energy state.

    Parameters
    ----------
    Z : int
        Atomic number.
    A : int
        Atomic mass.
    state : str
        energy state.

    Returns
    -------
    int
        Canonical nuclide id.

    Raises
    ------
    ValueError
        If the input energy state is invalid.

    Examples
    --------
    >>> rd.utils.build_id(1,2)
    10020000
    >>> rd.utils.build_id(28,56,'m')
    280560001

    """

    if state != "":
        if state in get_metastable_chars():
            state_int = get_metastable_chars().index(state) + 1
        else:
            raise ValueError(f"{state} is not a valid energy state.")
    else:
        state_int = 0

    canonical_id = (Z * 10000000) + (A * 10000) + state_int

    return canonical_id


def build_nuclide_string(Z: int, A: int, meta_state: str = "") -> str:
    """
    Builds a nuclide string from given atomic mass and number.

    Parameters
    ----------
    Z : int
        Atomic number.
    A : int
        Atomic mass.
    meta_state : str, optional
        Metastable state indicator character, ex. 'm' for the first
        atomic metastable state.

    Returns
    -------
    str
        Nuclide string built in symbol - mass number format.

    Raises
    ------
    ValueError
        If the input atomic number is invalid.

    Examples
    --------
    >>> rd.utils.build_nuclide_string(26, 56)
    'Fe-56'
    >>> rd.utils.build_nuclide_string(28, 56, 'm')
    'Fe-56m'

    """

    if Z not in Z_DICT:
        raise ValueError(f"{Z} is not a valid atomic number")

    return_string = f"{Z_DICT[Z]}-{A}{meta_state}"

    return return_string


class NuclideStrError(ValueError):
    """
    Custom exception class for invalid nuclide strings.

    Parameters
    ----------
    nuclide : str
        Nuclide string.
    additional_message : str
        Message with additional error context.

    """

    def __init__(self, nuclide: str, additional_message: str, *args) -> None:
        super().__init__(args)
        self.nuclide = nuclide
        self.additional_message = additional_message

    def __str__(self) -> str:
        return (
            f"{self.nuclide} is not a valid nuclide string. {self.additional_message}"
        )


def _process_metastable_element_str(metastable_element_str: str) -> Tuple[str, str]:
    """
    Function to separate a string that could either be an element string, or a combined metastable
    state char + element string.

    It separates the metastable state char correctly if the element is 2 chars (He, Li, Be, etc.).

    If the element is 1 char (H, B, C, N, O, F, I, U etc.), ambiguities can occur, e.g. is 'ni'
    Nickel or a second metastable state of Iodine. The atomic mass number & nuclear data is
    needed to correctly separate these cases. However for simplicity it is assumed user has
    capitalized correctly.

    Note: bugs may occur with incorrect capitalizations. See test_utils.py for some failure cases.
    """

    # metastable_element_str is 3 chars: must contain metastable state + 2 char element string
    if len(metastable_element_str) > 2:
        return metastable_element_str[0], metastable_element_str[1:]

    # metastable_element_str is 1 or 2 chars: assume metastable if first char is lower case and a
    # valid metastable state char, second char is uppercase and a valid element symbol
    if (
        metastable_element_str[0] in get_metastable_chars()
        and metastable_element_str[1:] in SYM_DICT
    ):
        return metastable_element_str[0], metastable_element_str[1:]

    return "", metastable_element_str  # Is ground state


def parse_nuclide_str(nuclide: str) -> str:
    """
    Parses a nuclide string from e.g. '241Pu' or 'Pu241' format to 'Pu-241' format. Note this
    function works for both radioactive and stable nuclides.

    Parameters
    ----------
    nuclide : str
        Nuclide string.

    Returns
    -------
    str
        Nuclide string parsed in symbol - mass number format.

    Raises
    ------
    NuclideStrError
        If the input nuclide string is invalid.

    Examples
    --------
    >>> rd.utils.parse_nuclide_str('222Rn')
    'Rn-222'
    >>> rd.utils.parse_nuclide_str('Ca40')
    'Ca-40'

    """

    original_input = nuclide
    nuclide = "".join(nuclide.split())  # Remove all whitespaces (Issue #65).
    nuclide = nuclide.replace("-", "", 1)  # Strip out first hyphen
    if not nuclide.isalnum():
        raise NuclideStrError(
            original_input,
            "Nuclide strings must contain only letters, numbers, and (optionally) up to one "
            "hyphen.",
        )

    A = "".join([n for n in nuclide if n.isdigit()])  # Mass number
    if len(A) == 0 or int(A) > 300:  # Largest mass number in NUBASE2020 is 295
        raise NuclideStrError(original_input, f"Mass number ({A}) is unphysical.")

    alpha_components = nuclide.split(A)
    if len(alpha_components) != 2:
        raise NuclideStrError(original_input, "")

    if alpha_components[0] == "":  # User inputted mass number first
        metastable_char, element = _process_metastable_element_str(alpha_components[1])
    else:  # User inputted element symbol first
        element = alpha_components[0]
        metastable_char = alpha_components[1]

    element = element.capitalize()
    if element not in SYM_DICT:
        raise NuclideStrError(original_input, f"Element ({element}) appears invalid.")

    if len(metastable_char) > 1:
        raise NuclideStrError(
            original_input,
            f"Ground state / metastable state specification ({metastable_char}) appears invalid.",
        )
    metastable_char = metastable_char.lower()
    if not (metastable_char == "" or metastable_char in get_metastable_chars()):
        raise NuclideStrError(
            original_input,
            f"Ground state / metastable state specification ({metastable_char}) appears invalid.",
        )

    return f"{element}-{A}{metastable_char}"


def parse_id(input_id: int) -> str:
    """
    Parses a nuclide canonical id in zzzaaammmm format into symbol -
    mass number format.

    Parameters
    ----------
    input_id : int
        Nuclide in canonical id format.

    Returns
    -------
    str
        Nuclide string parsed in symbol - mass number format.

    Examples
    --------
    >>> rd.utils.parse_id(190400000)
    'K-40'
    >>> rd.utils.parse_id(280560001)
    'Ni-56m'

    """

    id_zzzaaa = int(input_id / 10000)
    state_digits = input_id - (id_zzzaaa * 10000)
    state = get_metastable_chars()[state_digits - 1] if state_digits > 0 else ""
    Z = int(id_zzzaaa / 1000)
    A = id_zzzaaa - (Z * 1000)
    name = build_nuclide_string(Z, A, state)

    return name


def parse_nuclide(
    input_nuclide: Union[str, int], nuclides: np.ndarray, dataset_name: str
) -> str:
    """
    Parses a nuclide string or canonical id into symbol - mass number
    format and checks whether the nuclide is contained in the decay
    dataset.

    Parameters
    ----------
    input_nuclide : str or int
        Nuclide name string or canonical id in zzzaaammmm format.
    nuclides : numpy.ndarray
        List of all the nuclides in the decay dataset.
    dataset_name : str
        Name of the decay dataset.

    Returns
    -------
    str
        Nuclide string parsed in symbol - mass number format.

    Raises
    ------
    ValueError
        If the input nuclide string or id is invalid or the nuclide is
        not contained in the decay dataset.
    TypeError
        If the input is an invalid type, a string or integer is
        expected.

    Examples
    --------
    >>> rd.utils.parse_nuclide('222Rn', rd.DEFAULTDATA.nuclides, rd.DEFAULTDATA.dataset_name)
    'Rn-222'
    >>> rd.utils.parse_nuclide('Ba137m', rd.DEFAULTDATA.nuclides, rd.DEFAULTDATA.dataset_name)
    'Ba-137m'
    >>> rd.utils.parse_nuclide(280560001, rd.DEFAULTDATA.nuclides, rd.DEFAULTDATA.dataset_name)
    'Ni-56m'

    """

    original_input = input_nuclide
    if isinstance(input_nuclide, int):
        name = parse_id(input_nuclide)
    elif isinstance(input_nuclide, str):
        name = input_nuclide
    else:
        raise TypeError("Invalid input type, expected int or str")
    nuclide = parse_nuclide_str(name)

    if nuclide not in nuclides:
        raise ValueError(
            f"{original_input} is not a valid nuclide in {dataset_name} decay dataset."
        )

    return nuclide


def add_dictionaries(
    dict_a: Dict[str, Union[float, Expr]], dict_b: Dict[str, Union[float, Expr]]
) -> Dict[str, Union[float, Expr]]:
    """
    Adds together two dictionaries of nuclide strings and associated quantities. Supports both
    floats or SymPy quantities.

    Parameters
    ----------
    dict_a : dict
        First dictionary containing nuclide strings as keys and quantitites as values.
    dict_b : dict
        Second dictionary containing nuclide strings as keys and quantitites as values.

    Returns
    -------
    dict
        Combined dictionary containing the nuclides in both dict_a and dict_b, where the
        quantities have been added together when a nuclide is present in both input
        dictionaries.

    Examples
    --------
    >>> dict_a = {'Pm-141': 1.0, 'Rb-78': 2.0}
    >>> dict_b = {'Pm-141': 3.0, 'Rb-90': 4.0}
    >>> rd.utils.add_dictionaries(dict_a, dict_b)
    {'Pm-141': 4.0, 'Rb-78': 2.0, 'Rb-90': 4.0}

    """

    new_dict = dict_a.copy()
    for nuclide, quantity in dict_b.items():
        if nuclide in new_dict:
            new_dict[nuclide] = new_dict[nuclide] + quantity
        else:
            new_dict[nuclide] = quantity

    return new_dict


def sort_dictionary_alphabetically(
    input_inv_dict: Dict[str, Union[float, Expr]]
) -> Dict[str, Union[float, Expr]]:
    """
    Sorts a dictionary alphabetically by its keys.

    Parameters
    ----------
    input_inv_dict : dict
        Dictionary containing nuclide strings or Nuclide objects as keys and numbers
        as values.

    Returns
    -------
    dict
        Inventory dictionary which has been sorted by the nuclides alphabetically.

    Examples
    --------
    >>> rd.utils.sort_dictionary_alphabetically({'U-235': 1.2, 'Tc-99m': 2.3, 'Tc-99': 5.8})
    {'Tc-99': 5.8, 'Tc-99m': 2.3, 'U-235': 1.2}

    """

    return dict(sorted(input_inv_dict.items(), key=lambda x: x[0]))


def sort_list_according_to_dataset(
    input_list: List[str], key_dict: Dict[str, int]
) -> List[str]:
    """
    Sorts a list of nuclides based on their order of appearence in the decay dataset.

    Parameters
    ----------
    input_list : list
        List of nuclide strings to be sorted.
    key_dict : dict
        Dictionary from the decay dataset with nuclide strings as keys and their position
        (integers) in the decay dataset.

    Returns
    -------
    list
        Sorted nuclide list.

    Examples
    --------
    >>> rd.utils.sort_list_according_to_dataset(['Tc-99', 'Tc-99m'], rd.DEFAULTDATA.nuclide_dict)
    ['Tc-99m', 'Tc-99']

    """

    return sorted(input_list, key=lambda nuclide: key_dict[nuclide])
