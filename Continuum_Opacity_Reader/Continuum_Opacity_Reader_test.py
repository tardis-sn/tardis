"""This is the Continuum Opacity Source Reader test file. It will be used as a test file for the

Continuum Opacity Source Reader project where there will be an attempt to read from scientific literature and turn it to a

table that can be used by the project as a whole for TARDIS simulations.

Date : 2024-06-11

Author : Sanhik Nandi

Version : Prototype 0.01

"""

import pandas as pd
import gzip

class OpacityReader:
    """
    This is the reader class inside the Continuum Opacity Source Reader class. This class is currently used to take .gz
    files and return a Pandas dataframe that the TARDIS project can use to ingest data. It is basically a data processing
    pipeline for the project.
    The class has 2 methods:
    1. __init__(self, filepath, skip_rows, masking_value): this is the default constructor for the file and is used to
    initialize the variables used in the class. Returns nothing.
    2. read(): This is the core reader for the Continuum Opacity Source Reader project. This is the method used for
    the actual reading process and is used to work on the data. Returns the Pandas Dataframe.
    """
    def __init__(self, file_path, masking_value = None, skip_rows = 0):
        """
        Initializes the class Continuum Opacity Source Reader Class object with file and parsing parameters.
        Args:
            file_path (str): The path to the file to be read. Currently, the module only takes .gz files.

            masking_value(int, optional): The value used for masking the row values. A scientific paper has special rows that have values
            that describe the other rows instead of showing values. This parameter is used to mask those rows. By default,
            the masking value is set to None. Ensuring users get full access to class data.

            skip_rows (int, optional): The number of rows to skip from the beginning of the file, used to omit empty rows in the
            file which have little value and can skew the dataset. Defaults to 0 if no value is given to have the user see
            the whole dataset without any filtering.

        Returns: None

        """
        self.file_path = file_path
        self.masking_value = masking_value
        self.skip_rows = skip_rows

    def read(self):
        """Reads the data from the specified file path and returns it as a pandas DataFrame.

        Returns :
                Pandas Dataframe, containing the datatable extrapolated by the method.
        """
        try:
            column_names = ['index', 'log_t', 'log_r', 'opacity']
            with gzip.open(self.file_path, 'rt') as f:
                data = pd.read_csv(
                    f,
                    sep = r'\s+',
                    engine = 'python',
                    header = None,
                    skiprows = self.skip_rows
                )
            data.dropna(axis = 1, how = 'all', inplace = True)
            data.columns = column_names
            if (self.masking_value is not None):
                data['index']=pd.to_numeric(data['index'], errors = 'coerce')
                data = data.dropna(subset = ['index'])
                data = data[data['index'] < self.masking_value]

            return data
        except Exception as e:
            print(f"An error occurred while reading the file: {e}")
            return None

