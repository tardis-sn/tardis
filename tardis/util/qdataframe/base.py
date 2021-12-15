from typing import Type
import pandas as pd
from pandas.api.extensions import (
    ExtensionArray,
    ExtensionDtype,
    register_extension_dtype,
)

from astropy.units.quantity import Quantity
from collections.abc import Iterable 
import numpy as np

@register_extension_dtype
class QuantityDtype(ExtensionDtype):
    type=np.dtype
    name="quantity"
    unit=Quantity
    
    @classmethod
    def construct_array_type(cls):
        return QuantityArray
    
    @classmethod
    def construct_dtype_from(cls, quantity):
        dtype = cls()
        # dtype.type = quantity.dtype
        dtype.unit = quantity.unit
        return dtype
    


class QuantityArray(ExtensionArray):
    
    
    def __init__(self, data):
        self.data = data
        if not isinstance(self.data, Quantity):
            raise TypeError
        if not isinstance(self.data.value, Iterable):
            raise TypeError
        self._dtype = QuantityDtype.construct_dtype_from(self.data)
        
        
        
    @classmethod
    def _from_sequence(cls, scalars, dtype=None, copy=False, units=None):
        if isinstance(scalars, Quantity):
            if scalars.isscalar:
                return [scalars.value] * scalars.unit
            else:
                return scalars
        if not isinstance(scalars, Quantity) and units:
            # if a sequence and units are provided
            # return a sequence of quantities
            scalars = scalars * units
            return scalars
        else:
            raise ValueError
    

    @classmethod
    def _from_factorized(cls, values, original):
        return values * original.unit
    
    def __getitem__(self, idx):
        return self.data[idx]
    
    #TODO: Implement __setitem__
    def __setitem__(self, key, value):
        self.data.__setitem__(key, value)
    
    def __len__(self):
        return len(self.data)
    
    def __eq__(self, other):
        return self.data.__eq__(other)

    @property
    def dtype(self):
        return self._dtype
    
    @property
    def nbytes(self):
        return self._dtype.dtype.itemsize * len(self)
    
    def isna(self):
        return np.isnan(self.data)
    
    def take(self, indices, allow_fill=False, fill_value=None):
        from pandas.core.algorithms import take

        # If the ExtensionArray is backed by an ndarray, then
        # just pass that here instead of coercing to object.
        data = self.data

        if allow_fill and fill_value is None:
            fill_value = np.nan

        # fill value should always be translated from the scalar
        # type for the array, to the physical storage type for
        # the data, before passing to take.

        result = take(data, indices, fill_value=fill_value,
                    allow_fill=allow_fill)
        return self._from_sequence(result, dtype=self.dtype)
    
        
    def copy(self):
        return QuantityArray(self.data)
    
    @classmethod
    def _concat_same_type(cls, to_concat):
        """
        Concatenate multiple array

        Parameters
        ----------
        to_concat : sequence of this type

        Returns
        -------
        ExtensionArray
        """
        data = np.concatenate([ga.data for ga in to_concat])
        return QuantityArray(data)
    

class QuantityArray2(ExtensionArray, Quantity):
    def __new__(cls, *args, **kwargs):
        print("argsa", args)
        print("kwargs", kwargs)
        
        # return Quantity.__new__(Quantity, *args, **kwargs)
        res = super().__new__(Quantity, *args, **kwargs)
        print(res)
        return res
    
    def __init__(self, data):
        print("here")
        Quantity.__init__(self, data)
        print("here2")
        self.data = data
        if not isinstance(self.data, Quantity):
            raise TypeError
        if not isinstance(self.data.value, Iterable):
            raise TypeError
        self._dtype = QuantityDtype.construct_dtype_from(self.data)
        # self.dtype = self._dtype
        
        
        
    @classmethod
    def _from_sequence(cls, scalars, dtype=None, copy=False, units=None):
        if isinstance(scalars, Quantity):
            if scalars.isscalar:
                return [scalars.value] * scalars.unit
            else:
                return scalars
        if not isinstance(scalars, Quantity) and units:
            # if a sequence and units are provided
            # return a sequence of quantities
            scalars = scalars * units
            return scalars
        else:
            raise ValueError
    

    @classmethod
    def _from_factorized(cls, values, original):
        return values * original.unit
    
    def __getitem__(self, idx):
        return self.data[idx]
    
    #TODO: Implement __setitem__
    def __setitem__(self, key, value):
        self.data.__setitem__(key, value)
    
    def __len__(self):
        return len(self.data)
    
    def __eq__(self, other):
        return self.data.__eq__(other)

    @property
    def dtype(self):
        return self._dtype
    
    @property
    def nbytes(self):
        return self._dtype.dtype.itemsize * len(self)
    
    def isna(self):
        return np.isnan(self.data)
    
    def take(self, indices, allow_fill=False, fill_value=None):
        from pandas.core.algorithms import take

        # If the ExtensionArray is backed by an ndarray, then
        # just pass that here instead of coercing to object.
        data = self.data

        if allow_fill and fill_value is None:
            fill_value = np.nan

        # fill value should always be translated from the scalar
        # type for the array, to the physical storage type for
        # the data, before passing to take.

        result = take(data, indices, fill_value=fill_value,
                    allow_fill=allow_fill)
        return self._from_sequence(result, dtype=self.dtype)
    
        
    def copy(self):
        return QuantityArray(self.data)
    
    @classmethod
    def _concat_same_type(cls, to_concat):
        """
        Concatenate multiple array

        Parameters
        ----------
        to_concat : sequence of this type

        Returns
        -------
        ExtensionArray
        """
        data = np.concatenate([ga.data for ga in to_concat])
        return QuantityArray(data)
    
    def __repr__(self):
        from pandas.io.formats.printing import format_object_summary

        # the short repr has no trailing newline, while the truncated
        # repr does. So we include a newline in our template, and strip
        # any trailing newlines from format_object_summary
        # data = format_object_summary(
        #     self, self._formatter(), indent_for_name=False
        # ).rstrip(", \n")
        # class_name = f"<{type(self).__name__}>\n"
        # return f"{class_name}{data}\nLength: {len(self)}, dtype: {self.dtype}"
        return ExtensionArray.__repr__(self)
    
    def _formatter(self, boxed: bool = False):
        if boxed:
            return str
        return repr