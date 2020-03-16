import numpy as np
import pandas as pd
from pandas import DataFrame, Series
from astropy.units import Unit

from tardis.container.units import Units
from tardis.container.base import UnitContainer

def _validate_columns(df, units):
    if not units:
        return

    if isinstance(units, (dict, Series, Units)):
        key_list = units.keys()
        if units and not set(key_list).issubset(set(df.columns)):
            raise KeyError(f"Keys of {units} are not in DataFrame columns")
    else:
        raise TypeError(f"units cannot be {type(units)}")

class UDataFrame(DataFrame, UnitContainer):
    """
    Extended version of `Pandas.DataFrame`

    This data structure contains unit in each columns of `Pandas.DataFrame`.
    Each columns of unit cannot be exist.

    Parameters
    ----------
    data: `Numpy.ndarray`, Iterable, dict, or DataFrame
    index: `Pandas.Index` or array-like
    columns: `Pandas.Index` or array-like
    dtype: `Pandas.dtype`, default None
    copy: bool, default False
    units: array-like, Iterable, dict, or scalar value

    Notes
    ----------
    Most operations follow `Pandas.DataFrame`
    """

    _metadata = ["_units"]

    def __init__(self, *args, **kwargs):
        units = kwargs.pop("units", None)
        super(UDataFrame, self).__init__(*args, **kwargs)
        _validate_columns(self, units)
        self._units = Units(units, index=self.columns)

    @property
    def _constructor(self):
        return UDataFrame
    
    def copy(self, deep=True):
        """
        Make a copy of UDataFrame

        Units in UDataFrame also copy too.

        Parameters
        ----------
        deep : bool, default True

        Returns
        -------
        UDataFrame

        See Also
        --------
        More information needed, See Pandas Reference.

        DataFrame.copy : Make a copy of DataFrame
        """
        return super(UDataFrame, self).copy(deep=deep).__finalize__(self)

    def drop(self, *args, **kwargs):
        """
        Remove specific labels from columns or rows

        If specific labels are in columns, remove specific labels together.

        Parameters
        ----------
        labels : single label or list-like
        axis : {0 or 'index', 1 or 'columns'}, default 0
        index : single label or list-like
        columns : single label or list-like
        level : int or level name, optional
        inplace : bool, default False
        errors : {'ignore', 'raise'}, default 'raise'

        Returns
        -------
        UDataFrame

        Notes
        ------
        Drop Units together

        See Also
        --------
        More information needed, See Pandas Reference.

        DataFrame.drop : Remove specific label from columns or rows
        """
        inplace = kwargs.get("inplace", False)
        if inplace:
            super(UDataFrame, self).drop(*args, **kwargs)

            # Drop specific index in Units
            columns = kwargs.pop("columns", None)
            axis = kwargs.pop("axis", 0)
            if axis == 1:
                self.units.drop(*args, **kwargs)
            elif columns is not None:
                kwargs.update({
                    "index": columns
                })
                self.units.drop(*args, **kwargs)
            
            return self
        else:
            df = super(UDataFrame, self).drop(*args, **kwargs)

            # Drop specific index in Units
            columns = kwargs.pop("columns", None)
            axis = kwargs.pop("axis", 0)
            if axis == 1:
                _units = self.units.drop(*args, **kwargs)
            elif columns is not None:
                kwargs.update({
                    "index": columns
                })
                _units = self.units.drop(*args, **kwargs)
        
            df.units = _units if not _units.empty else self.units
            return df 
    
    def reindex(self, *args, **kwargs):
        """
        Change UDataFrame to new index with optional filling logic.

        Index of Units in UDataFrame also changes to new index 

        Parameters
        ----------
        labels : array-like, optional
        index, columns : array-like, optional
        axis: int or str, optional
        method : {None, 'backfill'/'bfill', 'pad'/'ffill', 'nearest'}
        copy : bool, default True
        level : int or name
        fill_value : scalar, default np.NaN
        limit : int, default None
        tolerance : optional

        Returns
        -------
        UDataFrame

        See Also
        --------
        More information needed, See Pandas Reference.

        DataFrame.reindex : Change UDataFrame to new index with optional filling logic.
        """
        df = super(UDataFrame, self).reindex(*args, **kwargs).__finalize__(self)

        # reindex Units
        axis = kwargs.pop("axis", None)
        columns = kwargs.pop("columns", None)
        copy = kwargs.pop("copy", True)
        level = kwargs.pop("level", None)
        if axis == 1 or axis == "columns":
            _units = self.units.reindex(*args, copy=copy, level=level)
        elif columns is not None:
            _units = self.units.reindex(columns, copy=copy, level=level)

        df.units = df.units if _units.empty else _units
        return df

    def rename(self, *args, **kwargs):
        """
        Change axes labels.

        Rename Units in UDataFrame together.

        Parameters
        ----------
        mapper : dict-like or function
        index : dict-like or function
        columns : dict-like or function
        axis : int or str
        copy : bool, default True
        inplace : bool, default False
        level : int or level name, default None
        errors : {'ignore', 'raise'}, default 'ignore'

        Returns
        -------
        UDataFrame

        Notes
        ------
        Rename Units together.

        See Also
        ---------
        More information needed, See Pandas Reference.

        DataFrame.rename: Change axes labels.
        """
        inplace = kwargs.get("inplace", None)
        if inplace:
            super(UDataFrame, self).rename(*args, **kwargs)

            # rename index in Units
            columns = kwargs.pop("columns", None)
            axis = kwargs.pop("axis", 0)
            if axis == 1:
                self.units.rename(*args, **kwargs)
            elif columns is not None:
                self.units.rename(columns, **kwargs)

            return self
        else:
            df = super(UDataFrame, self).rename(*args, **kwargs).__finalize__(self)

            # rename index in Units
            columns = kwargs.pop("columns", None)
            axis = kwargs.pop("axis", 0)
            if axis == 1:
                _units = df.units.rename(*args, **kwargs)
            elif columns is not None:
                _units = df.units.rename(columns, **kwargs)

            df.units = df.units if _units.empty else _units
            return df
    