"""
The fileio module contains functions to create inventory instances by reading nuclide quantities
from a file.

The docstring code examples assume that ``radioactivedecay`` has been imported
as ``rd``:

.. highlight:: python
.. code-block:: python

    >>> import radioactivedecay as rd

"""

import csv
import pathlib
from typing import Any, Dict, List, Optional, Tuple, Union

from radioactivedecay.decaydata import DecayData
from radioactivedecay.inventory import Inventory, InventoryHP


def _read_csv_file(
    filepath: Union[str, pathlib.Path],
    delimiter: str,
    encoding: str,
) -> List[List[str]]:
    """Read CSV file. All file read side-effect is here to assist testing read_csv() with mock."""
    with open(filepath, "r", encoding=encoding) as file:
        reader_object = csv.reader(file, delimiter=delimiter)
        lines = list(reader_object)

    return lines


def _parse_row(
    row: List[str], default_unit: Optional[str]
) -> Tuple[Dict[Union[str, int], float], Optional[str]]:
    """Parse one row that should be nuclide, quantity (, unit) format."""

    if len(row) not in {2, 3}:
        raise ValueError(
            f"This row of input file ('{row}') does not satisfy required format: i.e. nuclide, "
            "quantity (, unit)."
        )

    nuc = row[0]
    nuclide: Union[str, int] = int(nuc) if nuc.isnumeric() else nuc

    quantity = float(row[1])

    unit = default_unit
    if len(row) == 3:
        unit = row[2]

    return {nuclide: quantity}, unit


def read_csv(
    filepath: Union[str, pathlib.Path],
    inventory_type: str = "Inventory",
    units: Optional[str] = None,
    decay_data: Optional[DecayData] = None,
    delimiter: str = ",",
    skip_rows: int = 0,
    encoding: str = "utf-8",
) -> Union[Inventory, InventoryHP]:
    """
    Create an inventory by reading from a CSV file.

    The CSV file must contain at least two columns. The first column must contain nuclide strings /
    canonical ids. The second column must contain the quantities of each nuclide.

    Optionally a third column may be included specifying the units of each nuclide's quantity. If a
    unit is present on a row, it overrides the units set in the function parameter (i.e. ``units``
    is ignored for that row).

    A header row is not required. If the CSV file contains a header, you should ignore it using the
    ``skip_rows`` parameter.

    Example valid CSV file (no units - 'Bq' will be inferred by inventory constructor default):

    .. code-block:: text

        H-3,1.0
        C-14,2.0

    Example with tab separators (TSV format) and units specified:

    .. code-block:: text

        H-3 3.0 Bq
        He-3    17.0   mol
        C-14    5.0 Ci

    Parameters
    ----------
    filepath : str or pathlib.Path
        Name or path of the file to read in.
    inventory_type : str, optional
        The type of inventory you want to create. Either 'Inventory' (default) for a normal
        precision (SciPy) inventory, or 'InventoryHP' for a high precision (SymPy) inventory.
    units : None or str, optional
        The units of all the nuclide quantities in the file. If a unit is specified on a row of the
        file in the third column, that unit will override this one. If ``units = None`` is
        specified (default) and there is no unit on the row, the default of the inventory
        constructor will be used (currently 'Bq').
    decay_data : None or DecayData, optional
        The decay dataset to create the inventory with. If None is specified (default), the default
        of the inventory constructor will be used (which currently is the ICRP-107 dataset).
    delimiter : str, optional
        The delimiter used to separate items in the file. Default is comma (',', CSV format).
        Ex. if ``delimiter = '\t'`` (tab) will attempt to read a TSV file.
    skip_rows : int, optional
        Number of rows to ignore at the start of the file. Default is 0. Use this parameter
        to ignore a header row if the file has one (e.g. ``nuclide,quantity,unit``).
    encoding : str, optional
        Encoding of the file. Default is 'utf-8'.

    Raises
    ------
    ValueError
        If the file or function parameters are invalid in some way: e.g. i) if a row does not
        contain 2 (nuclide string/canonical id & quantity) or 3 values (if including unit), or ii)
        skip_rows means end of file is reached, or iii) invalid ``inventory_type`` specified.

    """

    lines = _read_csv_file(filepath, delimiter, encoding)
    lines_no_header = lines[skip_rows:]
    if len(lines_no_header) == 0:
        raise ValueError(
            f"`skip_rows = {skip_rows}`, but file only contains {len(lines)} lines so there cannot"
            " be any nuclide and quantity entries (rows) in the file."
        )

    contents, row_unit = _parse_row(lines_no_header[0], units)
    kwargs: Dict[str, Any] = {"contents": contents}
    if row_unit is not None or units is not None:
        kwargs["units"] = row_unit or units
    if decay_data is not None:
        kwargs["decay_data"] = decay_data

    if inventory_type == "Inventory":
        inv: Union[Inventory, InventoryHP] = Inventory(**kwargs)
    elif inventory_type == "InventoryHP":
        inv = InventoryHP(**kwargs)
    else:
        raise ValueError(
            f"`inventory_type` argument ('{inventory_type}') should either 'Inventory' or "
            "'InventoryHP'"
        )

    for row in lines_no_header[1:]:
        contents, row_unit = _parse_row(row, units)
        kwargs = {"add_contents": contents}
        if row_unit is not None or units is not None:
            kwargs["units"] = row_unit or units
        inv.add(**kwargs)

    return inv
