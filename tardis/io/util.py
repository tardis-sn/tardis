#Utility functions for the IO part of TARDIS

import pandas as pd
import numpy as np
import collections
from collections import OrderedDict
import yaml
import copy
from astropy import constants, units as u
from tardis.util import element_symbol2atomic_number

import logging
logger = logging.getLogger(__name__)


def quantity_from_str(text):
    """
    Convert a string to `astropy.units.Quantity`
    Parameters
    ----------
    text:
        The string to convert to `astropy.units.Quantity`
    Returns
    -------
    `astropy.units.Quantity`
    """
    value_str, unit = text.split(None, 1)
    value = float(value_str)
    if unit.strip() == 'log_lsun':
        value = 10 ** (value + np.log10(constants.L_sun.cgs.value))
        unit = 'erg/s'
    return u.Quantity(value, unit)


class MockRegexPattern(object):
    """
    A mock class to be used in place of a compiled regular expression
    when a type check is needed instead of a regex match.

    Note: This is usually a lot slower than regex matching.
    """
    def __init__(self, target_type):
        self.type = target_type

    def match(self, text):
        """

        Parameters
        ----------
        text:
            A string to be passed to `target_type` for conversion.
        Returns
        -------
        `True` if `text` can be converted to `target_type`.
        """
        try:
            self.type(text)
        except ValueError:
            return False
        return True


class YAMLLoader(yaml.Loader):
    """
    A custom YAML loader containing all the constructors required
    to properly parse the tardis configuration.
    """
    @classmethod
    def add_implicit_resolver(cls, tag, regexp, first=None):
        """
        Parameters
        ----------

        tag:
            The YAML tag to implicitly apply to any YAML scalar that matches `regexp`

        regexp:
            The regular expression to match YAML scalars against `tag`


        Notes
        -----

        This classmethod is a monkey-patch for a copy() related bug
        in the original class method which affects this yaml.Loader subclass.

        This class method is to be removed when this bug gets fixed upstream.

        https://bitbucket.org/xi/pyyaml/issues/57/add_implicit_resolver-on-a-subclass-may
        """
        if 'yaml_implicit_resolvers' not in cls.__dict__:
            yaml_implicit_resolvers = {}
            for k, v in cls.yaml_implicit_resolvers.items():
                yaml_implicit_resolvers[k] = copy.copy(v)
            cls.yaml_implicit_resolvers = yaml_implicit_resolvers
        if first is None:
            first = [None]
        for ch in first:
            cls.yaml_implicit_resolvers.setdefault(ch, []).append((tag, regexp))

    def construct_quantity(self, node):
        """
        A constructor for converting quantity-like YAML nodes to
        `astropy.units.Quantity` objects.

        Parameters
        ----------

        node:
            The YAML node to be constructed

        Returns
        -------

        `astropy.units.Quantity`

        """
        data = self.construct_scalar(node)
        return quantity_from_str(data)

    def mapping_constructor(self, node):
        return OrderedDict(self.construct_pairs(node))

YAMLLoader.add_constructor(u'!quantity', YAMLLoader.construct_quantity)
YAMLLoader.add_implicit_resolver(u'!quantity',
                                 MockRegexPattern(quantity_from_str))
YAMLLoader.add_implicit_resolver(u'tag:yaml.org,2002:float',
                                 MockRegexPattern(float))
YAMLLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
                           YAMLLoader.mapping_constructor)


def parse_abundance_dict_to_dataframe(abundance_dict):
    atomic_number_dict = dict([(element_symbol2atomic_number(symbol), abundance_dict[symbol])
                                   for symbol in abundance_dict])
    atomic_numbers = sorted(atomic_number_dict.keys())

    abundances = pd.Series([atomic_number_dict[z] for z in atomic_numbers], index=atomic_numbers)

    abundance_norm = abundances.sum()
    if abs(abundance_norm - 1) > 1e-12:
        logger.warn('Given abundances don\'t add up to 1 (value = %g) - normalizing', abundance_norm)
        abundances /= abundance_norm

    return abundances


def traverse_configs(base, other, func, *args):
    """
    Recursively traverse a base dict or list along with another one
    calling `func` for leafs of both objects.

    Parameters
    ----------
    base:
        The object on which the traversing is done
    other:
        The object which is traversed along with `base`
    func:
        A function called for each leaf of `base` and the correspnding leaf of `other`
        Signature: `func(item1, item2, *args)`
    args:
        Arguments passed into `func`

    """
    if isinstance(base, collections.Mapping):
        for k in base:
            traverse_configs(base[k], other[k], func, *args)
    elif isinstance(base, collections.Iterable) and not isinstance(base, basestring) and not hasattr(base, 'shape'):
        for val1, val2 in zip(base, other):
            traverse_configs(val1, val2, func, *args)
    else:
        func(base, other, *args)


def assert_equality(item1, item2):
    assert type(item1) is type(item2)
    try:
        if hasattr(item1, 'unit'):
            assert item1.unit == item2.unit
        assert np.allclose(item1, item2, atol=0.0)
    except (ValueError, TypeError):
        assert item1 == item2


def check_equality(item1, item2):
    try:
        traverse_configs(item1, item2, assert_equality)
    except AssertionError:
        return False
    else:
        return True
