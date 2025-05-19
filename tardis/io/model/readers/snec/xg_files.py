from io import StringIO
from pathlib import Path
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.units import Quantity
import xarray as xr

@dataclass
class XGData:
    timestamps: Quantity
    data_blocks: list
    metadata: dict = field(default_factory=dict)

    def to_xarray(self):
        """
        Converts the XGData to an xarray DataArray.

        Returns
        -------
        xr.DataArray: A 3D xarray DataArray where each DataFrame in data_blocks is a slice along the first dimension.
        """

        # Ensure all DataFrames share the same index and columns
        idx = self.data_blocks[0].index
        cols = self.data_blocks[0].columns

        # Stack their values into a single numpy array
        arr = np.stack([df.values for df in self.data_blocks], axis=0)

        # Create the xarray DataArray
        da3d = xr.DataArray(
            arr,
            coords={
                'time': self.timestamps,
                'cell_id': idx+1,
                'quantity': cols
            },
            dims=('time', 'cell_id', 'quantity')
        )

        return da3d



def xg_block_size(path):
    """
    Determines the block size (number of lines) between consecutive lines starting with \"Time in an .xg file.

    Parameters
    ----------
        path (str): Path to the .xg file.

    Returns
    -------
        int: The block size.
    """
    timestamp_blocks = []

    with open(path) as fh:
        for i, line in enumerate(fh):
            if line.startswith(' "Time'):
                timestamp_blocks.append(i)

    # Ensure all differences between timestamps are consistent
    assert (
        len(timestamp_blocks) > 1
    ), "File must contain at least two 'Time' entries to determine block size."
    assert np.all(
        np.diff(timestamp_blocks) == np.diff(timestamp_blocks)[0]
    ), "Inconsistent block sizes detected in the file."

    return len(timestamp_blocks), np.diff(timestamp_blocks)[0]


def read_xg_file(
    file_path: str, column_names: list, show_progress: bool = False
):
    """
    Reads the timestamps and corresponding data blocks from an .xg file.

    Parameters
    ----------
        file_path (str): Path to the .xg file.
        column_names (list): List of column names for the data blocks.
        show_progress (bool): Whether to show a progress bar (default: False).

    Returns
    -------
        list: A list of timestamps extracted from the file.
        list: A list of pandas DataFrames, each containing a data block corresponding to a timestamp.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"The file {file_path} does not exist.")

    timestamps = []
    data_blocks = []

    # Read the entire file into memory
    with open(path) as file_obj:
        lines = file_obj.readlines()

    # Split the file into blocks based on lines starting with "Time
    num_blocks, block_size = xg_block_size(path)

    block_ranges = list(
        zip(
            np.arange(0, num_blocks) * block_size,
            np.arange(1, num_blocks + 1) * block_size,
        )
    )

    iterator = block_ranges
    if show_progress:
        try:
            from IPython.core import get_ipython

            ip = get_ipython()
            if ip is not None and "IPKernelApp" in ip.config:
                from tqdm.notebook import tqdm as tqdm_progress
            else:
                from tqdm import tqdm as tqdm_progress
        except ImportError:
            from tqdm import tqdm as tqdm_progress
        iterator = tqdm_progress(iterator, total=num_blocks, desc="Reading xg blocks")

    for block_start, block_end in iterator:
        # Extract the timestamp from the current block
        current_block = lines[block_start:block_end]
        timestamp = float(current_block[0].split("=")[1].strip().split()[0])
        timestamps.append(timestamp)

        # Extract the data block corresponding to the timestamp
        data_block = pd.read_csv(
            StringIO("".join(current_block[1:])), sep=r"\s+", header=None,
            names=column_names
        )
        data_blocks.append(data_block)

    return XGData(timestamps * u.s, data_blocks)
