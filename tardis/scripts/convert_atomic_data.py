"""
Quick-and-dirty script for converting an older atomic data format to one compatible with modern TARDIS.
"""

import argparse
import hashlib
import pickle
import platform
import uuid
from datetime import datetime

import h5py
import numpy as np
import pandas as pd
import pytz


def serialize_pandas_object(pd_object):
    """
    Serialize Pandas objects with Pickle.

    Parameters
    ----------
    pd_object : pandas.Series or pandas.DataFrame
        Pandas object to be serialized with Pickle.

    Returns
    -------
    Pickle serialized Python object.
    """
    return pickle.dumps(pd_object)


def hash_pandas_object(pd_object, algorithm="md5"):
    """
    Hash Pandas objects.

    Parameters
    ----------
    pd_object : pandas.Series or pandas.DataFrame
        Pandas object to be hashed.
    algorithm : str, optional
        Algorithm available in `hashlib`, by default "md5"

    Returns
    -------
    str
        Hash values.

    Raises
    ------
    ValueError
        If `algorithm` is not available in `hashlib`.
    """
    algorithm = algorithm.lower()

    if hasattr(hashlib, algorithm):
        hash_func = getattr(hashlib, algorithm)

    else:
        raise ValueError('algorithm not supported')

    return hash_func(serialize_pandas_object(pd_object)).hexdigest()


def multiindex_port(olddata, templatedata, templatekey, oldkey=None):
    if oldkey is None:
        oldkey = templatekey

    # Get the format of the multi-index from the desired template
    index_names = templatedata[templatekey].index.names

    # Attempt conversion of the old structured array to a pd DataFrame
    newdata = pd.DataFrame(olddata[oldkey][:]).set_index(index_names)

    # Check datatypes of columns, convert if necessary
    if hasattr(templatedata[templatekey], "columns"):
        for col in templatedata[templatekey].columns:
            desired_dtype = templatedata[templatekey][col].dtype
            if desired_dtype == np.dtype("object"):
                desired_dtype = str
            newdata[col] = newdata[col].astype(desired_dtype)
    else:
        col = [templatedata[templatekey].name]
        newdata[col] = newdata[col].astype(templatedata[templatekey].dtype)

    # Convert datatypes of each level of the multi-index
    if not isinstance(newdata.index, pd.MultiIndex):
        template_ind_dtype = templatedata[templatekey].index.dtype
        newdata.index = newdata.index.astype(template_ind_dtype)
    else:
        for i, indname in enumerate(newdata.index.names):
            template_ind_dtype = templatedata[templatekey].index.dtypes[indname]
            converted_index = newdata.index.levels[i].astype(template_ind_dtype)
            newdata.index = newdata.index.set_levels(converted_index, level=i)
    return newdata


def simple_port(olddata, templatedata, templatekey, oldkey=None):
    if oldkey is None:
        oldkey = templatekey

    newdata = pd.DataFrame(olddata[oldkey][:])

    # Check datatypes of columns, convert if necessary
    for col in templatedata[templatekey].columns:
        newdata[col] = newdata[col].astype(templatedata[templatekey][col].dtype)

    return newdata


## Files used to convert Christian's atomic data:
#default_oldatomdata = "merged_mod_20SNG_forbidden_yg_fix_H30_cmfgen_yg.h5"
#default_pi_filename = "photoionization_data_H30_He.h5"
#default_template = "/home/connor/tardis-regression-data/atom_data/nlte_atom_data/TestNLTE_He_Ti.h5"
#default_new = "test2.h5"

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--oldatomdata", type=str)#, default=default_oldatomdata)
    parser.add_argument("--oldPIdata", type=str)#, default=default_pi_filename)
    parser.add_argument("--template", type=str)#, default=default_template)
    parser.add_argument("--newatomdata", type=str)#, default=default_new)
    args = parser.parse_args()

    # Atomic data in old format - "/Y/g" seems to be the only accessible dataset when
    # using pandas for some reason (contains collision data)
    old_df = pd.HDFStore(args.oldatomdata)

    # Photoionization data stored in separate file in Christian's version
    pi_data = pd.HDFStore(args.oldPIdata)

    # Reference atomic data file for the format we want to convert to
    template = pd.HDFStore(args.template)

    # Open up a new pandas HDFStore to port the old data into
    new = pd.HDFStore(args.newatomdata)

    ### COLLISIONS DATA
    multiindex_cols = list(old_df["/Y/g"].columns[:4])
    new["collisions_data"] = old_df["/Y/g"].set_index(multiindex_cols)
    tempcols = list(new['collisions_data'].columns)
    new["collisions_data_temperatures"] = pd.Series(tempcols, dtype=np.int64)
    new['collisions_data'] = new['collisions_data'].rename(lambda f: tempcols.index(f), axis=1)

    # We've gotten all useful info out of the old dataframe using pandas
    # so now load the same file with h5py directly to get the rest
    old_df.close()
    old = h5py.File(args.oldatomdata)

    ### COLLISIONS METADATA
    new_metadata = template["collisions_metadata"].copy()
    new_metadata["temperatures"] = new['collisions_data_temperatures'].values
    new["collisions_metadata"] = new_metadata

    ### ATOM DATA
    new_atom_data = multiindex_port(
        old, template, "atom_data", oldkey="basic_atom_data"
    )
    new["atom_data"] = new_atom_data

    ### IONIZATION DATA
    new["ionization_data"] = multiindex_port(old, template, "ionization_data")

    ### LEVELS DATA
    new["levels_data"] = multiindex_port(old, template, "levels_data")

    ### LINES DATA
    new["lines_data"] = multiindex_port(old, template, "lines_data")

    ### LINES METADATA
    new["lines_metadata"] = template["lines_metadata"].copy()

    ### MACRO ATOM DATA
    new["macro_atom_data"] = simple_port(old, template, "macro_atom_data")

    ### MACRO ATOM REFERENCES
    new["macro_atom_references"] = multiindex_port(
        old, template, "macro_atom_references"
    )

    ### PHOTOIONIZATION DATA
    ### Note this comes from a DIFFERENT file!
    new["photoionization_data"] = pi_data["photoionization_data"].copy()

    ### ZETA DATA
    zeta_index = pd.MultiIndex.from_arrays(
        old["zeta_data"][:, :2].T.astype(np.int64),
        names=template["zeta_data"].index.names,
    )
    zeta_temps = pd.Index(
        old["zeta_data"].attrs["t_rad"].astype(np.float64), name="temp"
    )
    new_zeta_data = pd.DataFrame(
        old["zeta_data"][:, 2:], index=zeta_index, columns=zeta_temps
    )
    new["zeta_data"] = new_zeta_data

    ### METADATA
    # Copied over from Andrew's notebook demonstrating how to do this
    meta = []
    meta.append(("format", "version", "2.0"))

    total_checksum = hashlib.md5()
    for key in new.keys():
        # update the total checksum to sign the file
        total_checksum.update(serialize_pandas_object(new[key]))

        # save individual DataFrame/Series checksum
        checksum = hash_pandas_object(new[key])
        meta.append(("md5sum", key.lstrip("/"), checksum))

    # relevant package versions
    meta.append(("software", "python", platform.python_version()))
    imports = [
        "carsus",
        "astropy",
        "numpy",
        "pandas",
        "tables",
        "ChiantiPy",
    ]
    for package in imports:
        meta.append(("software", package, __import__(package).__version__))
    meta_df = pd.DataFrame.from_records(
        meta, columns=["field", "key", "value"], index=["field", "key"]
    )
    uuid1 = uuid.uuid1().hex
    new.root._v_attrs["MD5"] = total_checksum.hexdigest()
    new.root._v_attrs["UUID1"] = uuid1
    new.root._v_attrs["FORMAT_VERSION"] = "2.0"
    tz = pytz.timezone("UTC")
    date = datetime.now(tz).isoformat()
    new.root._v_attrs["DATE"] = date
    new.put("/metadata", meta_df)

    # Close / write out files
    old.close()
    template.close()
    pi_data.close()
    new.close()
