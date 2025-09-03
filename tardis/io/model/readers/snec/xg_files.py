from dataclasses import dataclass, field
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from astropy import units as u
from astropy.units import Quantity
from tqdm.auto import tqdm


@dataclass
class XGData:
    timestamps: Quantity
    data_blocks: list
    metadata: dict = field(default_factory=dict)

    def to_xr_dataset(self):
        idx = self.data_blocks[0].index
        cols = [
            col
            for col in self.data_blocks[0].columns
            if col not in ["enclosed_mass", "radius"]
        ]
        arr = np.stack(
            [
                df.drop(["enclosed_mass", "radius"], axis=1).values
                for df in self.data_blocks
            ],
            axis=0,
        )

        # Enclosed mass stays 1D
        enclosed_mass = self.data_blocks[0]["enclosed_mass"].values

        # Radius becomes 2D (time, cell_id)
        radius_arr = np.stack([df["radius"].values for df in self.data_blocks], axis=0)

        data_vars = {
            col: (("time", "cell_id"), arr[:, :, i]) for i, col in enumerate(cols)
        }

        ds = xr.Dataset(
            data_vars=data_vars,
            coords={"time": self.timestamps, "cell_id": idx + 1},
        )
        ds = ds.assign_coords(enclosed_mass=("cell_id", enclosed_mass))
        ds = ds.assign_coords(radius=(("time", "cell_id"), radius_arr))

        return ds


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


def read_xg_file(file_path: str, column_names: list, show_progress: bool = True):
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

    iterator = tqdm(block_ranges) if show_progress else block_ranges
    for block_start, block_end in iterator:
        # Extract the timestamp from the current block
        current_block = lines[block_start:block_end]
        timestamp = float(current_block[0].split("=")[1].strip().split()[0])
        timestamps.append(timestamp)

        # Extract the data block corresponding to the timestamp
        data_block = pd.read_csv(
            StringIO("".join(current_block[1:])),
            sep=r"\s+",
            header=None,
            names=column_names,
        )
        data_blocks.append(data_block)

    return XGData(timestamps * u.s, data_blocks)
