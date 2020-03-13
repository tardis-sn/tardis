import pandas as pd
import re
import os

"""
Python program to read the tabulated data present in si2_osc_kurucz file,
parse them using regular expressions, and tabulate into dataframes.
"""

class ParameterException(Exception):
    pass

def find(name, path):
    """
    Helper function to iteratively find the required file in the directory.

    Parameters
    ----------
    name : name of file.

    path : absolute path of the directory.

    Returns
    -------
    Absolute path of the file.
    """
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)

def extract_data(path):
    """
    Reads the data present in the file and return as a string.

    Parameters
    ----------
    path : Absolute path to the file.
    
    Returns
    -------
    Filedata (in string format).
    """
    filedata = ""
    with open(path, 'r') as f:
        filedata += f.read()
    return filedata

def extract_tabular_data(pattern, string, columns, flags=0):
    """
    Extracts the tabulated data from the file according to the pattern.

    Parameters
    ----------
    pattern : Regex pattern for matching.

    string : string on which to match the expression.

    columns : list of column names.

    flags : (Optional) flags to be passed to matcher object.

    Returns
    -------
    Pandas Dataframe containing the data.
    """
    df = pd.DataFrame(columns=columns)
    matches = re.finditer(pattern, string, flags)

    for matchNum, match in enumerate(matches, start=1):
        for groupNum in range(0, len(match.groups())):
            df.loc[matchNum, columns[groupNum]] = match.group(groupNum+1)
    
    return df

def read_data(directory_path, filename, pattern, columns):
    """
    Driver function.
    """
    # Finds file and extracts data
    file_path = find(filename, directory_path)
    if file_path is None:
        raise ParameterException('File not found at speified location.')
    print('File found at path:', file_path)

    filedata = extract_data(file_path)
    
    # Metadata parsing
    df = extract_tabular_data(pattern, filedata, columns, re.MULTILINE)
    return df
    # Extraction of tabulated data
    df1 = extract_tabular_data(PATTERN1, table1, columns1, re.MULTILINE)
    df1.set_index("ID", inplace=True)
    print(df1.head())

    df2 = extract_tabular_data(PATTERN2, table2, columns2, re.MULTILINE)
    print(df2.head())