import pandas as pd
from pandas import DataFrame, Series
from astropy.units import Unit

from tardis.container.units import Units

def _delegate_method_to_units(op, this, *args, **kwargs):
    return getattr(this, op)(*args, **kwargs)

class UnitContainer(object):

    @property
    def units(self):
        """
        Units represented as a extended version of`Pandas.Series`

        Returns
        -------
        Units
        """
        return self._units
    
    @units.setter
    def units(self, value):
        """
        Set units

        Notes
        -----
        Keys of value must be equal to current keys of units
        """
        value = value if isinstance(value, Units) else Units(value)
        self._units = value
    
    def get_unit(self, label):
        """
        Get unit of single column.

        Parameters
        -----------
        label : column label

        Returns
        --------
        Unit
        """
        return _delegate_method_to_units("get_unit", self.units, label=label)
    
    def set_unit(self, label, unit):
        """
        Put unit at passed column

        Parameters
        ----------
        label : column label

        Returns
        -------
        Unit
        """
        return _delegate_method_to_units("set_unit", self.units, label=label, unit=unit)

    def copy(self, *args, **kwargs):
        """
        Make a copy of object.
        
        See Also
        --------
        See Pandas Documentation 

        DataFrame.copy : Make a copy of DataFrame
        Series.copy : Make a copy of Series
        """
        raise NotImplementedError()

    def drop(self, *args, **kwargs):
        """
        Drop specific labels from rows or columns.

        If drop UDataFrame columns, then drop columns in Units together.

        See also
        ---------
        See Pandas Documentation

        DataFrame.drop : Return DataFrame with specific index or column labels removed.
        Series.drop : Return Series with specific index labels removed.
        """
        raise NotImplementedError()

    def reindex(self, *args, **kwargs):
        """
        Change UnitContainer to new index with optional filling logic.

        If reindex UDataFrame columns, then reindex Units together.

        See also
        --------
        See Pandas Documentation

        DataFrame.reindex : Change DataFrame to new index with optional filling logic.
        Series.reindex : Change Series to new index with optional filling logic.
        """
        raise NotImplementedError()

    def rename(self, *args, **kwargs):
        """
        Change axes labels.

        See Also
        ---------
        DataFrame.rename : Return DataFrame with specified index 
            or column labels renamed.
        Series.rename : Return Series with specified index labels renamed.
        """
        raise NotImplementedError()

