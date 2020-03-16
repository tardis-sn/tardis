import numpy as np
import pandas as pd
from pandas import DataFrame, Series, isna
from pandas.core.dtypes.common import is_scalar
from astropy.units import Unit


class Units(Series):
    """
    Container of list of Units

    This data structure contains astropy.units objects which can be accessed by index.
    Most operation same as `Pandas.Series`.
    """

    def __init__(self, *args, **kwargs):
        if not kwargs and not args:
            raise ValueError("Cannot construct empty Units")

        super(Units, self).__init__(*args, **kwargs)

        for k, v in self.items():
            if is_scalar(v) and isna(v):
                self.loc[k] = v
            else:
                try:
                    self.loc[k] = Unit(v)
                except (TypeError, ValueError):
                    raise TypeError(f"Cannot type cast to Unit")

    @property
    def _constructor(self):
        return Units

    @property
    def units(self):
        """
        Units represented as a `Pandas.Series`

        Returns
        -------
        Series

        Notes
        -----
        If no units, empty Series
        """
        return self
    
    @units.setter
    def units(self, value):
        """
        Set units from value

        Notes
        -----
        Keys of value must be equal to current keys of units
        """
        value = value if isinstance(value, Units) else Units(value)
        self = value        

    def get_unit(self, label):
        """
        Get Unit from specific label

        Parameters
        ----------
        label : Single label or list-like

        Returns
        -------
        Units
        """
        if label not in self.index:
            raise KeyError(f"{label} not in index")
        else:
            return self.loc[label]
    
    def set_unit(self, label, unit):
        """
        Set Unit to specific label

        Parameters
        ----------
        label : Single label or list-like
        unit : astropy.units or str

        Returns
        -------
        Units
        """
        if label not in self.index:
            raise KeyError(f"{label} not in index")
        else:
            try:
                u = Unit(unit)
            except (TypeError, ValueError):
                raise TypeError(f"{unit} cannot be cast to astropy.units.Unit")

            self.loc[label] = u

    def copy(self, deep=True):
        """
        Make a copy of Units

        Parameters
        ----------
        deep : bool, default True

        Returns
        -------
        Units

        See Also
        --------
        More information needed, See Pandas Reference.

        Series.copy : Make a copy of Series
        """
        return super(Units, self).copy(deep=deep)

    def drop(self, *args, **kwargs):
        """
        Remove specific labels from rows.

        Parameters
        ----------
        labels : single label or list-like
        axis : 0, default 0
        index : single label or list-like
        columns : single label or list-like
        level : int or level name, optional
        inplace : bool, default False
        errors : {'ignore', 'raise'}, default 'raise'

        Returns
        -------
        Units

        See Also
        ---------
        More information needed, See Pandas Reference.

        Series.drop: Remove specific label from rows
        """
        inplace = kwargs.get("inplace", False)
        if inplace:
            super(Units, self).drop(*args, **kwargs)
            return self
        else:
            return super(Units, self).drop(*args, **kwargs)

    def reindex(self, index, copy=True, level=None, fill_value=np.NaN):
        """
        Change Units to new index

        This parameters of operation are different from Series.reindex.
        This operation doesn't offer optional filling logic, 
        only can fill with scala value

        Parameters
        ----------
        index : array-like, optional
        copy : bool, default True
        level : int or name
        fill_value : scalar, default np.NaN

        Returns
        -------
        Units

        See Also
        --------
        More information needed, See Pandas Reference.

        Series.reindex : Change Series to new index with optional filling logic.
        """
        return super(Units, self).reindex(
            index, copy=copy, level=level, fill_value=fill_value
        )

    def rename(self, *args, **kwargs):
        """
        Change axes labels

        Parameters
        ----------
        axis : {0 or "index"}
        index : scalar, hashable sequence, dict-like or function, optional
        inplace : bool, default False

        Returns
        -------
        Units
        
        See Also
        ---------
        More information needed, See Pandas Reference.

        Series.rename: Change index labels.
        """
        inplace = kwargs.get("inplace", False)
        if inplace:
            super(Units, self).rename(*args, **kwargs)
            return self
        else:
            return super(Units, self).rename(*args, **kwargs)
