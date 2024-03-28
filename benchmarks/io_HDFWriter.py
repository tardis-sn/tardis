"""
Basic TARDIS Benchmark.
"""
import numpy as np
import pandas as pd
from astropy import units as u
from asv_runner.benchmarks.mark import parameterize, skip_benchmark

from benchmarks.benchmark_base import BenchmarkBase
from tardis.io.util import HDFWriterMixin


# @skip_benchmark
class BenchmarkIoHdfWriter(BenchmarkBase):
    """
    Class to benchmark the HDF Writer function.
    """

    def __init__(self):
        arr = np.array(
            [
                ["L1", "L1", "L2", "L2", "L3", "L3", "L4", "L4"],
                ["one", "two", "one", "two", "one", "two", "one", "two"],
            ]
        )
        self.mock_multiIndex = pd.MultiIndex.from_arrays(arr.transpose())

    class MockHDF(HDFWriterMixin):
        hdf_properties = ["property"]

        def __init__(self, property):
            self.property = property

    simple_objects = [1.5, "random_string", 4.2e7]

    @parameterize({"Simple objects": simple_objects})
    def time_simple_write(self, attr):
        fname = self.hdf_file_path
        actual = self.MockHDF(attr)
        actual.to_hdf(fname, path="test", overwrite=True)
        pd.read_hdf(fname, key="/test/mock_hdf/scalars")

    mock_df = pd.DataFrame(
        {
            "one": pd.Series([1.0, 2.0, 3.0], index=["a", "b", "c"]),
            "two": pd.Series([1.0, 2.0, 3.0, 4.0], index=["a", "b", "c", "d"]),
        }
    )

    complex_objects = [
        np.array([4.0e14, 2, 2e14, 27.5]),
        pd.Series([1.0, 2.0, 3.0]),
        mock_df,
    ]

    @parameterize({"Complex objects": complex_objects})
    def time_complex_obj_write(self, attr):
        fname = self.hdf_file_path
        actual = self.MockHDF(attr)
        actual.to_hdf(fname, path="test", overwrite=True)
        pd.read_hdf(fname, key="/test/mock_hdf/property")

    def time_multi_index_write(self):
        fname = self.hdf_file_path
        actual = self.MockHDF(self.mock_multiIndex)
        actual.to_hdf(fname, path="test", overwrite=True)
        expected = pd.read_hdf(fname, key="/test/mock_hdf/property")
        pd.MultiIndex.from_tuples(expected.unstack().values)

    # Test Quantity Objects
    quantity_objects = [np.array([4.0e14, 2, 2e14, 27.5]), mock_df]

    @parameterize({"Quantity objects": quantity_objects})
    def time_quantity_objects_write(self, attr):
        fname = self.hdf_file_path
        attr_quantity = u.Quantity(attr, "g/cm**3")
        actual = self.MockHDF(attr_quantity)
        actual.to_hdf(fname, path="test", overwrite=True)
        pd.read_hdf(fname, key="/test/mock_hdf/property")

    scalar_quantity_objects = [1.5, 4.2e7]

    @parameterize({"Scalar quantity objects": scalar_quantity_objects})
    def time_scalar_quantity_objects_write(self, attr):
        fname = self.hdf_file_path
        attr_quantity = u.Quantity(attr, "g/cm**3")
        actual = self.MockHDF(attr_quantity)
        actual.to_hdf(fname, path="test", overwrite=True)
        pd.read_hdf(fname, key="/test/mock_hdf/scalars/")

    def time_none_write(self):
        fname = self.hdf_file_path
        actual = self.MockHDF(None)
        actual.to_hdf(fname, path="test", overwrite=True)
        pd.read_hdf(fname, key="/test/mock_hdf/scalars/")

    # Test class_properties parameter (like homologous_density is a class
    # instance/object inside Model class)
    class MockClass(HDFWriterMixin):
        hdf_properties = ["property", "nested_object"]

        def __init__(self, property, nested_object):
            self.property = property
            self.nested_object = nested_object

    @parameterize({"Quantity objects": quantity_objects})
    def time_objects_write(self, attr):
        fname = self.hdf_file_path
        nested_object = self.MockHDF(np.array([4.0e14, 2, 2e14, 27.5]))
        attr_quantity = u.Quantity(attr, "g/cm**3")
        actual = self.MockClass(attr_quantity, nested_object)
        actual.to_hdf(fname, path="test", overwrite=True)
        pd.read_hdf(fname, key="/test/mock_class/property")
        pd.read_hdf(
            fname, key="/test/mock_class/nested_object/property"
        )

    def time_snake_case(self):
        assert (
                self.MockHDF.convert_to_snake_case("HomologousDensity")
                == "homologous_density"
        )
        assert self.MockHDF.convert_to_snake_case("TARDISSpectrum") == "tardis_spectrum"
        assert self.MockHDF.convert_to_snake_case("BasePlasma") == "base_plasma"
        assert self.MockHDF.convert_to_snake_case("LTEPlasma") == "lte_plasma"
        assert (
                self.MockHDF.convert_to_snake_case("MontecarloTransport")
                == "montecarlo_transport"
        )
        assert (
                self.MockHDF.convert_to_snake_case("homologous_density")
                == "homologous_density"
        )
