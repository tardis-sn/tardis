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
    'tardis/visualization/__init__.py',
    'tardis/visualization/widgets/__init__.py',
]:
    write_package(_pkg, '')

write_package('tardis/visualization/__init__.py', """\
from tardis.visualization.widgets.shell_info import shell_info_from_simulation
import panel as pn
pn.extension("tabulator")
""")

# stub of tardis.util.base, which is needed by the widgets
write_package('tardis/util/base.py', """\
import pandas as pd
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

def atomic_number2element_symbol(atomic_number):
    return ATOMIC_NUMBER2SYMBOL[atomic_number]

def species_tuple_to_string(species_tuple, roman_numerals=True):
    atomic_number, ion_number = species_tuple
    element_symbol = ATOMIC_NUMBER2SYMBOL[atomic_number]
    if roman_numerals:
        roman_ion_number = int_to_roman(ion_number + 1)
        return f"{str(element_symbol)} {roman_ion_number}"
    return f"{element_symbol} {ion_number:d}"
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


resp = await pyodide.http.pyfetch(f'{GITHUB_URL}tardis/visualization/widgets/util.py')
write_package('tardis/visualization/widgets/util.py', await resp.string())

resp = await pyodide.http.pyfetch(f'{GITHUB_URL}tardis/visualization/widgets/shell_info.py')
write_package('tardis/visualization/widgets/shell_info.py', await resp.string())

resp = await pyodide.http.pyfetch(STATIC_URL + 'data/atomic_symbols.dat')
write_data('/data/atomic_symbols.dat', await resp.string())

resp = await pyodide.http.pyfetch(STATIC_URL + 'data/shell_info_data.pkl')
write_data('/data/shell_info_data.pkl', await resp.bytes())
