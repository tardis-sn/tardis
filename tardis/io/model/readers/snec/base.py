from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from tqdm import tqdm


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


def read_xg_timestamps_and_data(
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
            from IPython import get_ipython

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
            StringIO("".join(current_block[1:])), sep=r"\s+", header=None
        )
        data_blocks.append(data_block)

    return timestamps * u.s, data_blocks
