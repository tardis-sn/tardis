from astropy.units import Unit
import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal, assert_series_equal
import pytest

from tardis.container.uframe import UDataFrame
from tardis.container.units import Units

class TestUDataFrameConstructor:
    def test_constructor_with_unit(self):
        result = UDataFrame(
            {"A": [1, 2], "B": [3, 4], "C": [5, 6]},
            units={"A": "kg", "B": "meter", "C": "second"}
        )
        expected = UDataFrame({"A": [1, 2], "B": [3, 4], "C": [5, 6]})
        expected.units = Units(["kg", "meter", "second"], index=["A", "B", "C"])

        assert_frame_equal(result, expected)
        assert_series_equal(result.units, expected.units)

    def test_constructor_not_matching_columns_with_unit(self):
        msg = "are not in DataFrame columns"
        with pytest.raises(KeyError, match=msg):
            UDataFrame(
                {"A": [1, 2], "B": [3, 4], "C": [5, 6]},
                units={"A": "kg", "b": "meter", "D": "second"}
            )
    
    def test_constructor_empty_units(self):
        result = UDataFrame({"A": [1, 2], "B": [3, 4], "C": [5, 6]})
        expected = Units(None, index=["A", "B", "C"])

        assert_series_equal(result.units, expected)

class TestUDataFrameAttribute:
    def test_get_units(self):
        udf = UDataFrame(
            {"A": [1, 2], "B": [3, 4], "C": [5, 6]},
            units={"A": "kg", "B": "meter", "C": "second"}
        )

        result = udf.units
        expected = Units({"A": "kg", "B": "meter", "C": "second"})
        assert_series_equal(result, expected)

    def test_set_units(self):
        udf = UDataFrame(
            {"A": [1, 2], "B": [3, 4], "C": [5, 6]},
            units={"A": "kg", "B": "meter", "C": "second"}
        )

        expected = Units({"A": "Hz", "B": "mol", "C": "J"})
        udf.units = expected
        result = udf.units
        assert_series_equal(result, expected)

    def test_get_unit(self):
        udf = UDataFrame(
            {"A": [1, 2], "B": [3, 4], "C": [5, 6]},
            units={"A": "kg", "B": "meter", "C": "second"}
        )

        expected = Unit("meter")
        result = udf.get_unit("B")
        assert result == expected

        with pytest.raises(KeyError, match="not in index"):
            udf.get_unit("D")

    def test_set_unit(self):
        udf = UDataFrame(
            {"A": [1, 2], "B": [3, 4], "C": [5, 6]},
            units={"A": "kg", "B": "meter", "C": "second"}
        )

        expected = Units({"A": "Hz", "B": "mol", "C": "J"})
        udf.units = expected
        result = udf.units
        assert_series_equal(result, expected)

        with pytest.raises(KeyError, match="not in index"):
            udf.set_unit("D", expected)


class TestUDataFrameConversion:
    def test_copy(self):
        udf = UDataFrame(
            {"A": [1, 2], "B": [3, 4], "C": [5, 6]},
            units={"A": "kg", "B": "meter", "C": "second"}
        )

        result = udf.copy(deep=False)
        assert_frame_equal(udf, result)
        assert_series_equal(udf.units, result.units)

        result = udf.copy(deep=True)
        assert_frame_equal(udf, result)
        assert_series_equal(udf.units, result.units)

class TestUDataFrameReindexing:
    def test_drop_column(self):
        udf = UDataFrame(
            {"A": [1, 2], "B": [3, 4], "C": [5, 6]},
            units={"A": "kg", "B": "meter", "C": "second"}
        )
        expected = UDataFrame(
            {"B": [3, 4], "C": [5, 6]}, units={"B": "meter", "C": "second"}
        )
        
        result = udf.drop(columns=["A"])
        assert_frame_equal(result, expected)
        assert_series_equal(result.units, expected.units)

        result = udf.drop(["A"], axis=1)
        assert_frame_equal(result, expected)
        assert_series_equal(result.units, expected.units)

    def test_drop_column_with_inplace(self):
        udf = UDataFrame(
            {"A": [1, 2], "B": [3, 4], "C": [5, 6]},
            units={"A": "kg", "B": "meter", "C": "second"}
        )
        expected = UDataFrame({"C": [5, 6]}, units={"C": "second"})

        result = udf.drop(columns=["A", "B"], inplace=False)
        assert_frame_equal(result, expected)
        assert_series_equal(result.units, expected.units)        

        udf.drop(columns=["A", "B"], inplace=True)
        assert_frame_equal(udf, expected)
        assert_series_equal(udf.units, expected.units)

    def test_reindex_columns(self):
        udf = UDataFrame(
            {"A": [1, 2], "B": [3, 4], "C": [5, 6]},
            units={"A": "kg", "B": "meter", "C": "second"}
        )
        expected = UDataFrame(
            {"A": [1, 2], "B": [3, 4], "D": [np.NaN, np.NaN]},
            units={"A": "kg", "B": "meter", "D": np.NaN}
        )

        result = udf.reindex(columns=["A", "B", "D"])
        assert_frame_equal(result, expected)
        assert_series_equal(result.units, expected.units)

        result = udf.reindex(["A", "B", "D"], axis=1)
        assert_frame_equal(result, expected)
        assert_series_equal(result.units, expected.units)

    def test_reindex_columns_wiht_fill_value(self):
        udf = UDataFrame(
            {"A": [1, 2], "B": [3, 4], "C": [5, 6]},
            units={"A": "kg", "B": "meter", "C": "second"}
        )
        expected = UDataFrame(
            {"A": [1, 2], "B": [3, 4], "D": [10, 10]},
            units={"A": "kg", "B": "meter", "D": np.NaN}
        )

        result = udf.reindex(columns=["A", "B", "D"], fill_value=10)
        assert_frame_equal(result, expected)
        assert_series_equal(result.units, expected.units)

    def test_rename_columns(self):
        udf = UDataFrame(
            {"A": [1, 2], "B": [3, 4], "C": [5, 6]},
            units={"A": "kg", "B": "meter", "C": "second"}
        )
        expected = UDataFrame(
            {"D": [1, 2], "E": [3, 4], "F": [5, 6]},
            units={"D": "kg", "E": "meter", "F": "second"}
        )
        
        result = udf.rename(columns={"A": "D", "B": "E", "C": "F"})
        assert_frame_equal(result, expected)
        assert_series_equal(result.units, expected.units)

        result = udf.rename({"A": "D", "B": "E", "C": "F"}, axis=1)
        assert_frame_equal(result, expected)
        assert_series_equal(result.units, expected.units)

    def test_rename_columns_with_inplace(self):
        udf = UDataFrame(
            {"A": [1, 2], "B": [3, 4], "C": [5, 6]},
            units={"A": "kg", "B": "meter", "C": "second"}
        )
        expected = UDataFrame(
            {"a": [1, 2], "b": [3, 4], "c": [5, 6]},
            units={"a": "kg", "b": "meter", "c": "second"}
        )

        result = udf.rename(str.lower, axis=1, inplace=False)
        assert_frame_equal(result, expected)
        assert_series_equal(result.units, expected.units)

        udf.rename(str.lower, axis=1, inplace=True)
        assert_frame_equal(udf, expected)
        assert_series_equal(udf.units, expected.units)
