"""
Setup code to run when any pyodide cell is executed.
This is responsible for:
    1. Creating lightweight mock modules for certain tardis modules so that
        pyodide does not need to install the whole dependency tree of tardis.
    2. Writing the files needed by pyodide into the virtual filesystem.

Note:
    This file is executed in a Pyodide web worker using pyodide.runPythonAsync().
    During execution, modules such as pyodide.http and js are available.
    The code runs asynchronously.
"""

import sysconfig
import pyodide.http
import js
from pathlib import Path


lib = Path(sysconfig.get_path('purelib'))  # Path to site-packages in pyodide's virtual filesystem

# Temporarily, this points to my fork
# Can be changed to the main repo later
GITHUB_URL = 'https://raw.githubusercontent.com/SS-9098/tardis/refs/heads/First-Objective-Migrate-to-Panel/'
# Derive _static/ base URL from the WebWorker's own script URL.
# The worker is always loaded from _static/PyodideWebWorker.js.
WORKER_URL = js.self.location.href
STATIC_URL = WORKER_URL[: WORKER_URL.rfind('/') + 1]  # trailing slash

def write_package(rel_path, content):
    full = lib / rel_path
    full.parent.mkdir(parents=True, exist_ok=True)
    full.write_text(content)


def write_data(rel_path, content):
    p = Path(rel_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    if p.suffix == '.pkl':
        p.write_bytes(content)
    else:
        p.write_text(content)


# Package __init__.py files
for _pkg in [
    'tardis/__init__.py',
    'tardis/util/__init__.py',
    'tardis/configuration/__init__.py',
    'tardis/transport/__init__.py',
    'tardis/transport/montecarlo/__init__.py',
    'tardis/transport/montecarlo/packets/__init__.py',
    'tardis/visualization/__init__.py',
    'tardis/visualization/widgets/__init__.py',
]:
    write_package(_pkg, '')

write_package('tardis/visualization/__init__.py', """\
from tardis.visualization.widgets.shell_info import shell_info_from_simulation
from tardis.visualization.widgets.line_info import LineInfoWidget
import panel as pn
pn.extension("tabulator", "plotly")
""")

# stub of tardis.util.base, which is needed by the widgets
write_package('tardis/util/base.py', """\
import pandas as pd
import re
from collections import OrderedDict

ATOMIC_SYMBOLS_DATA = (
    pd.read_csv(
        "/data/atomic_symbols.dat",
        sep=r"\s+",
        names=["atomic_number", "symbol"],
    )
    .set_index("atomic_number")
    .squeeze()
)

ATOMIC_NUMBER2SYMBOL = OrderedDict(ATOMIC_SYMBOLS_DATA.to_dict())
SYMBOL2ATOMIC_NUMBER = OrderedDict((y, x) for x, y in ATOMIC_NUMBER2SYMBOL.items())

NUMERAL_MAP = tuple(
    zip(
        (1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1),
        ("M", "CM", "D", "CD", "C", "XC", "L", "XL", "X", "IX", "V", "IV", "I"),
    )
)

def int_to_roman(i):
    result = []
    for integer, numeral in NUMERAL_MAP:
        count = i // integer
        result.append(numeral * count)
        i -= integer * count
    return "".join(result)

def roman_to_int(roman_string):
    NUMERALS_SET = set(list(zip(*NUMERAL_MAP))[1])
    roman_string = roman_string.upper()
    if len(set(roman_string.upper()) - NUMERALS_SET) != 0:
        raise ValueError(f"{roman_string} does not seem to be a roman numeral")
    i = result = 0
    for integer, numeral in NUMERAL_MAP:
        while roman_string[i : i + len(numeral)] == numeral:
            result += integer
            i += len(numeral)
    if result < 1:
        raise ValueError(f"Can not interpret Roman Numeral {roman_string}")
    return result

def reformat_element_symbol(element_string):
    return element_string[0].upper() + element_string[1:].lower()

def atomic_number2element_symbol(atomic_number):
    return ATOMIC_NUMBER2SYMBOL[atomic_number]
    
def element_symbol2atomic_number(element_string):
    reformatted_element_string = reformat_element_symbol(element_string)
    return SYMBOL2ATOMIC_NUMBER[reformatted_element_string]

def species_tuple_to_string(species_tuple, roman_numerals=True):
    atomic_number, ion_number = species_tuple
    element_symbol = ATOMIC_NUMBER2SYMBOL[atomic_number]
    if roman_numerals:
        roman_ion_number = int_to_roman(ion_number + 1)
        return f"{str(element_symbol)} {roman_ion_number}"
    return f"{element_symbol} {ion_number:d}"
    
def species_string_to_tuple(species_string):
    try:
        element_symbol, ion_number_string = re.match(
            r"^(\w+)\s*(\d+)", species_string
        ).groups()
    except AttributeError:
        try:
            element_symbol, ion_number_string = species_string.split()
        except ValueError:
            print(
                f'Species string "{species_string}" is not of format <element_symbol><number>'
                f" (e.g. Fe 2, Fe2, ..)"
            )

    atomic_number = element_symbol2atomic_number(element_symbol)

    try:
        ion_number = roman_to_int(ion_number_string)
    except ValueError:
        print(
            "Ion Number does not contain a Roman Numeral. Checking for integer value"
        )
        try:
            ion_number = int(ion_number_string)
        except ValueError:
            print(
                f"Given ion number ('{ion_number_string}') could not be parsed"
            )

    if ion_number - 1 > atomic_number:
        print("Species given does not exist: ion number > atomic number")

    return atomic_number, ion_number - 1
""")

write_package('tardis/util/environment.py', """\
class Environment:
    @staticmethod
    def allows_widget_display():
        return True

    @classmethod
    def get_current_environment(cls):
        return 'pyodide'
""")

write_package('tardis/transport/montecarlo/packets/radiative_packet.py', """\
from enum import IntEnum
class InteractionType(IntEnum):
    NO_INTERACTION = -1
    BOUNDARY = 1
    LINE = 2
    ESCATTERING = 4
    CONTINUUM_PROCESS = 8
""")

write_package('tardis/configuration/sorting_globals.py', """\
SORTING_ALGORITHM = "stable"
""")

write_package('tardis/constants.py', """\
from astropy.constants.astropyconst13 import *
""")

resp = await pyodide.http.pyfetch(f'{GITHUB_URL}tardis/visualization/widgets/util.py')
write_package('tardis/visualization/widgets/util.py', await resp.string())

resp = await pyodide.http.pyfetch(f'{GITHUB_URL}tardis/visualization/widgets/shell_info.py')
write_package('tardis/visualization/widgets/shell_info.py', await resp.string())

resp = await pyodide.http.pyfetch(f'{GITHUB_URL}tardis/analysis.py')
write_package('tardis/analysis.py', await resp.string())

resp = await pyodide.http.pyfetch(f'{GITHUB_URL}tardis/visualization/widgets/line_info.py')
write_package('tardis/visualization/widgets/line_info.py', await resp.string())

resp = await pyodide.http.pyfetch(f'{GITHUB_URL}tardis/visualization/plot_util.py')
write_package('tardis/visualization/plot_util.py', await resp.string())

resp = await pyodide.http.pyfetch(STATIC_URL + 'data/atomic_symbols.dat')
write_data('/data/atomic_symbols.dat', await resp.string())

resp = await pyodide.http.pyfetch(STATIC_URL + 'data/line_info_data.pkl')
write_data('/data/line_info_data.pkl', await resp.bytes())

resp = await pyodide.http.pyfetch(STATIC_URL + 'data/shell_info_data.pkl')
write_data('/data/shell_info_data.pkl', await resp.bytes())
