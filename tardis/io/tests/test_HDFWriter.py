import os

import numpy as np
import pandas as pd
import pandas.util.testing as pdt
import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from numpy.testing import assert_almost_equal, assert_array_almost_equal

from tardis.io.util import HDFWriter


#Test Cases

#DataFrame
#None
#Numpy Arrays
#Strings
#Numeric Values
#Pandas Series Object
#MultiIndex Object
#Quantity Objects with - Numeric Values, Numpy Arrays, DataFrame, Pandas Series, None objects

class MockHDF(HDFWriter, object):
    hdf_properties = ['property']
    class_properties = {}

    def __init__(self, property):
        self.property = property

simple_objects = [1.5, 'random_string', 4.2e7]

@pytest.mark.parametrize("attr", simple_objects)
def test_simple_write(tmpdir, attr):
    fname = str(tmpdir.mkdir('data').join('test.hdf'))
    actual = MockHDF(attr)
    actual.to_hdf(fname, path='test')
    expected = pd.read_hdf(fname, key='/test/mock_hdf/scalars')['property']
    assert actual.property == expected

mock_df = pd.DataFrame({'one': pd.Series([1., 2., 3.], index=['a', 'b', 'c']),
                        'two': pd.Series([1., 2., 3., 4.], index=['a', 'b', 'c', 'd'])})
complex_objects = [np.array([4.0e14, 2, 2e14, 27.5]),
                   pd.Series([1., 2., 3.]), mock_df]

@pytest.mark.parametrize("attr", complex_objects)
def test_complex_obj_write(tmpdir, attr):
    fname = str(tmpdir.mkdir('data').join('test.hdf'))
    actual = MockHDF(attr)
    actual.to_hdf(fname, path='test')
    expected = pd.read_hdf(fname, key='/test/mock_hdf/property').values
    assert_array_almost_equal(actual.property, expected)

arrays = [['L1', 'L1', 'L2', 'L2', 'L3', 'L3', 'L4', 'L4'],
          ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']]
tuples = list(zip(*arrays))
mock_multiIndex = pd.MultiIndex.from_tuples(tuples, names=['first', 'second'])

def test_MultiIndex_write(tmpdir):
    fname = str(tmpdir.mkdir('data').join('test.hdf'))
    actual = MockHDF(mock_multiIndex)
    actual.to_hdf(fname, path='test')
    expected = pd.read_hdf(fname, key='/test/mock_hdf/property')
    expected = pd.MultiIndex.from_tuples(
        expected.unstack().values, names=['first', 'second'])
    pdt.assert_almost_equal(actual.property, expected)

#Test Quantity Objects

quantity_objects = [np.array([4.0e14, 2, 2e14, 27.5]), mock_df]

@pytest.mark.parametrize("attr", quantity_objects)
def test_quantity_objects_write(tmpdir, attr):
    fname = str(tmpdir.mkdir('data').join('test.hdf'))
    attr_quantity = u.Quantity(attr, 'g/cm**3')
    actual = MockHDF(attr_quantity)
    actual.to_hdf(fname, path='test')
    expected = pd.read_hdf(fname, key='/test/mock_hdf/property')
    assert_array_almost_equal(actual.property.cgs.value, expected)

scalar_quantity_objects = [1.5, 4.2e7]

@pytest.mark.parametrize("attr", scalar_quantity_objects)
def test_scalar_quantity_objects_write(tmpdir, attr):
    fname = str(tmpdir.mkdir('data').join('test.hdf'))
    attr_quantity = u.Quantity(attr, 'g/cm**3')
    actual = MockHDF(attr_quantity)
    actual.to_hdf(fname, path='test')
    expected = pd.read_hdf(fname, key='/test/mock_hdf/scalars/')['property']
    assert_array_almost_equal(actual.property.cgs.value, expected)

def test_none_write(tmpdir):
    fname = str(tmpdir.mkdir('data').join('test.hdf'))
    actual = MockHDF(None)
    actual.to_hdf(fname, path='test')
    expected = pd.read_hdf(fname, key='/test/mock_hdf/scalars/')['property']
    if expected == 'none':
        expected = None
    assert actual.property == expected

# Test class_properties parameter (like homologous_density is a class
# instance/object inside Model class)

class MockClass(HDFWriter, object):
    hdf_properties = ['property', 'nested_object']
    class_properties = {'nested_object': MockHDF}

    def __init__(self, property, nested_object):
        self.property = property
        self.nested_object = nested_object

@pytest.mark.parametrize("attr", quantity_objects)
def test_objects_write(tmpdir, attr):
    fname = str(tmpdir.mkdir('data').join('test.hdf'))
    nested_object = MockHDF(np.array([4.0e14, 2, 2e14, 27.5]))
    attr_quantity = u.Quantity(attr, 'g/cm**3')
    actual = MockClass(attr_quantity, nested_object)
    actual.to_hdf(fname, path='test')
    expected_property = pd.read_hdf(fname, key='/test/mock_class/property')
    assert_array_almost_equal(actual.property.cgs.value, expected_property)
    nested_property = pd.read_hdf(
        fname, key='/test/mock_class/nested_object/property')
    assert_array_almost_equal(
        actual.nested_object.property, nested_property)
