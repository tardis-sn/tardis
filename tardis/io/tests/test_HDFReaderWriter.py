import os

import numpy as np
import pandas as pd
import pandas.util.testing as pdt
import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from numpy.testing import assert_almost_equal, assert_array_almost_equal

from tardis.io.util import HDFReaderWriter


class MockHDF(HDFReaderWriter, object):
    hdf_properties = ['property']
    quantity_attrs = {}
    class_properties = {}

    def __init__(self, property):
        self.property = property

#Test Cases

#DataFrame
#None
#Numpy Arrays
#Strings
#Numeric Values
#Pandas Series Object
#MultiIndex Object
#Quantity Objects with - Numeric Values, Numpy Arrays, DataFrame, Pandas Series, None objects


simple_objects = [1.5, 'random_string', 4.2e7, None]


@pytest.mark.parametrize("attr", simple_objects)
def test_simple_readwrite(tmpdir, attr):
    fname = str(tmpdir.mkdir('data').join('test.hdf'))
    actual = MockHDF(attr)
    actual.to_hdf(fname, path='mockHDF')
    expected = MockHDF.from_hdf(fname, path='mockHDF')
    assert actual.property == expected.property


mock_df = pd.DataFrame({'one': pd.Series([1., 2., 3.], index=['a', 'b', 'c']),
                        'two': pd.Series([1., 2., 3., 4.], index=['a', 'b', 'c', 'd'])})
complex_objects = [np.array([4.0e14, 2, 2e14, 27.5]),
                   pd.Series([1., 2., 3.]), mock_df]


@pytest.mark.parametrize("attr", complex_objects)
def test_complex_obj_readwrite(tmpdir, attr):
    fname = str(tmpdir.mkdir('data').join('test.hdf'))
    actual = MockHDF(attr)
    actual.to_hdf(fname, path='mockHDF')
    expected = MockHDF.from_hdf(fname, path='mockHDF')
    assert_array_almost_equal(actual.property, expected.property)


arrays = [['L1', 'L1', 'L2', 'L2', 'L3', 'L3', 'L4', 'L4'],
          ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']]
tuples = list(zip(*arrays))
mock_multiIndex = pd.MultiIndex.from_tuples(tuples, names=['first', 'second'])


def test_MultiIndex_readwrite(tmpdir):
    fname = str(tmpdir.mkdir('data').join('test.hdf'))
    actual = MockHDF(mock_multiIndex)
    actual.to_hdf(fname, path='mockHDF')
    expected = MockHDF.from_hdf(fname, path='mockHDF')
    expected.property = pd.MultiIndex.from_tuples(
        expected.property.unstack().values, names=['first', 'second'])
    pdt.assert_almost_equal(actual.property, expected.property)


quantity_objects = [np.array([4.0e14, 2, 2e14, 27.5]), mock_df, 1.5, 4.2e7]
quantity_attrs = {'property': 'g/cm**3'}


@pytest.mark.parametrize("attr", quantity_objects)
def test_quantity_objects_readwrite(tmpdir, attr):
    fname = str(tmpdir.mkdir('data').join('test.hdf'))
    actual = MockHDF(attr)
    actual.quantity_attrs = quantity_attrs
    actual.to_hdf(fname, path='mockHDF')
    expected = MockHDF.from_hdf(fname, path='mockHDF')
    assert_quantity_allclose(actual.property, expected.property)


def test_none_with_quantity_readwrite(tmpdir):
    fname = str(tmpdir.mkdir('data').join('test.hdf'))
    actual = MockHDF(None)
    actual.quantity_attrs = quantity_attrs
    actual.to_hdf(fname, path='mockHDF')
    expected = MockHDF.from_hdf(fname, path='mockHDF')
    assert actual.property == expected.property


class MockClass(HDFReaderWriter, object):
    hdf_properties = ['property']
    quantity_attrs = {}
    class_properties = {'nested_object': MockHDF}

    def __init__(self, property, nested_object):
        self.property = property
        self.nested_object = nested_object


# Test class_properties parameter (like homologous_density is a class instance/object inside Model class)
quantity_objects = [np.array([4.0e14, 2, 2e14, 27.5]), mock_df, 1.5, 4.2e7]
quantity_attrs = {'property': 'g/cm**3'}


@pytest.mark.parametrize("attr", quantity_objects)
def test_objects_readwrite(tmpdir, attr):
    fname = str(tmpdir.mkdir('data').join('test.hdf'))
    nested_object = MockHDF(np.array([4.0e14, 2, 2e14, 27.5]))
    actual = MockClass(attr, nested_object)
    actual.quantity_attrs = quantity_attrs
    actual.to_hdf(fname, path='mockHDF')
    expected = MockClass.from_hdf(fname, path='mockHDF')
    assert_quantity_allclose(actual.property, expected.property)
    assert_quantity_allclose(
        actual.nested_object.property, expected.nested_object.property)
