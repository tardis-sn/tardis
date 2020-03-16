from astropy.units import Unit
import numpy as np
import pandas as pd
from pandas.testing import assert_series_equal
import pytest

from tardis.container.units import Units

class TestUnitsConstructor:
    def test_constructor(self):
        result = Units(["kg", "meter", "second"], index=["A", "B", "C"])
        expected = Units(
            [Unit("kg"), Unit("meter"), Unit("second")],
            index=["A", "B", "C"]
        )
        assert_series_equal(result, expected)

        msg = "Cannot type cast to Unit"
        with pytest.raises(TypeError, match=msg):
            Units(["kg", "meter", "foo"], index=["A", "B", "C"])
    
class TestUnitsAttribute:
    def test_get_unit(self):
        result = Units(["kg", "meter", "second"], index=["A", "B", "C"])
        assert result.get_unit("B") == Unit("meter")

        with pytest.raises(KeyError, match="not in index"):
            result.get_unit("D")

    def test_set_unit(self):
        units = Units(["kg", "meter", "second"], index=["A", "B", "C"])
        
        expected = Unit("mol")
        units.set_unit("B", expected)
        assert units.get_unit("B") == expected

        with pytest.raises(KeyError, match="not in index"):
            units.set_unit("D", Unit("Hz"))


class TestUnitsConversion:
    def test_copy(self):
        units = Units(["kg", "meter", "second"], index=["A", "B", "C"])

        result = units.copy(deep=False)
        assert_series_equal(units, result)

        result = units.copy(deep=True)
        assert_series_equal(units, result)


class TestUnitsReindexing:
    def test_drop(self):
        units = Units(["kg", "meter", "second"], index=["A", "B", "C"])
        expected = Units(["meter", "second"], index=["B", "C"])
        
        result = units.drop("A", inplace=False)
        assert_series_equal(result, expected)

        units.drop("A", inplace=True)
        assert_series_equal(result, expected)
        
    def test_reindex(self):
        units = Units(["kg", "meter", "second"], index=["A", "B", "C"])

        msg = r"reindex\(\) got an unexpected keyword argument"
        with pytest.raises(TypeError, match=msg):
            Units.reindex(["A", "C", "D"], method="ffill")

    def test_rename(self):
        units = Units(["kg", "meter", "second"], index=["A", "B", "C"])
        new_index = {"A": "a", "C": "c"}
        expected = Units(["kg", "meter", "second"], index=["a", "B", "c"])

        result = units.rename(new_index, inplace=False)
        assert_series_equal(result, expected)

        units.rename(new_index, inplace=True)
        assert_series_equal(result, expected)
