"""Helpers to write compressed regression data during tests.

Use these functions in tests that write temporary regression HDF/NPY data so the
artifacts are smaller on CI and when generating reference data.
"""
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def save_npz(path: Path, arr: Any) -> None:
    """Save numpy arrays (or list of arrays) compressed using np.savez_compressed.

    Parameters
    ----------
    path: Path
        Destination file path, will get a .npz suffix if missing.
    arr: array-like
        A numpy array or a sequence of arrays to save.
    """
    path = Path(path)
    if path.suffix != ".npz":
        path = path.with_suffix(path.suffix + ".npz")
    if isinstance(arr, (list, tuple)):
        np.savez_compressed(path, *arr)
    else:
        np.savez_compressed(path, arr=arr)


def to_hdf_compressed(
    hdf_file_path: Path,
    obj: Any,
    key: str | None = None,
    name: str | None = None,
    overwrite: bool = True,
    format: str | None = None,
    *,
    complevel: int = 9,
    complib: str = "zlib",
) -> None:
    """Write pandas-friendly objects to HDF5 using compression where possible.

    This wraps pandas.DataFrame/Series.to_hdf or pd.HDFStore to set the
    complevel and complib arguments when a table format is used.

    Parameters
    ----------
    hdf_file_path: Path
        Destination HDF5 path.
    obj: Any
        Object to store. If it has a `to_hdf` method it will be used, otherwise
        pandas.DataFrame/Series will be tried.
    key: str
        HDF key under which to store the object.
    kwargs: dict
        Passed to the underlying `to_hdf` call.
    """
    hdf_file_path = Path(hdf_file_path)
    hdf_file_path.parent.mkdir(parents=True, exist_ok=True)

    # If the object has a to_hdf method (many domain objects do), prefer it.
    # The domain mixin already applies compression (blosc). We avoid passing
    # pandas compression args here to keep compatibility with the mixin signature.
    if hasattr(obj, "to_hdf"):
        try:
            group_name = name or key
            obj.to_hdf(hdf_file_path, name=group_name, overwrite=overwrite, format=format)
            return
        except TypeError:
            # Some implementations may differ; fall back to pandas path below
            pass

    # Fall back to pandas functions
    try:
        pd_obj = pd.DataFrame(obj)
        store_key = name or key or "/"
        pd_obj.to_hdf(
            hdf_file_path,
            key=store_key,
            format="table",
            complevel=complevel,
            complib=complib,
        )
    except Exception:
        # Last resort: use HDFStore directly without compression kwargs
        store_key = name or key or "/"
        store = pd.HDFStore(hdf_file_path, mode="a")
        store.put(store_key, pd.DataFrame(obj))
        store.close()
