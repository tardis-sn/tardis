import numpy as np

from pandas import DataFrame, Series

from astropy.units import Unit

class UDataFrame(DataFrame):
    """
    An extension of the `~pandas.Dataframe`, which stores the units of
    the columns as well and preserve during basic operations.

    Parameters
    ----------
    data : The actual data to be stored in dataframe

    units : pandas.Series of either astropy.units or None
            It stores the units corresponding to the columns of the data
    
    Examples
    --------
    >>> from tardis.model.udata import Udataframe
    >>> from pandas import DataFrame, Series
    >>> udata = UDataFrame(data=DataFrame(data = pd.np.random.randint(0, 50, (10, 5)),
                           columns=['A','B','C','D','E'],
                           units=Series([u.meter, None, u.kg, None, u.Ohm])
    """

    _metadata = ['units']

    @property
    def _constructor(self):
        return UDataFrame

    def __init__(self, data, units, **kwargs):
        units = Series(kwargs.pop("units", None))
        super().__init__(**kwargs)
        self.units = units

    def _verify_len_unit(self):
        """
        Verify if all the columns have an entry in the
        Series of units.
        """
        return True if len(self.columns) == len(self.units) else False

    def _verify_type_units(self):
        """
        Verify if all elements of units are either of type
        `~astropy.units.Unit` or None
        """
        return True if all(type(self.units) is Unit or None) else False 

    @property
    def units(self):
        return self._units

    @units.setter
    def units(self, u):
        if self._verify_len_unit and self._verify_type_units:
            if type(u) is Series:
                self._units = u
            elif type(u) is list:
                self._units = Series(u)
            else:
                raise TypeError("Units must be either a pandas.Series or a list")
        else:
            raise ValueError("The units must be of same length as that of data,",
                             " and elements must be either of type astropy.units.Unit or None")

    def get_unit(self, column_name):
        """
        Function to get unit of a column

        Parameters
        ----------
        column_name : str
               Column Name
        
        Returns
        -------
        `~astropy.units.Unit` or None
        """
        if column_name not in self.columns:
            raise ValueError("Invalid column name")
        
        i = list(self.columns).index(column_name)
        return self.units[i]

    def set_unit(self, column_name, unit_name):
        """
        Function to set unit of a column

        Parameters
        ----------
        column_name : str
               Column Name
        unit_name : `~astropy.units.Unit` or None
                    Unit to be set to the column
        """
        if column_name not in self.columns:
            raise ValueError("Invalid column name")
        
        i = list(self.columns).index(column_name)
        self.units[i] = unit_name

    def copy(self):
        """
        Creates a copy of the original dataframe

        Returns
        -------
        `~tardis.model.udata.UDataFrame
        """
        cp = UDataFrame(self)
        cp.units = self.units
        return cp
