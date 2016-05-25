#Utility functions for the IO part of TARDIS

import pandas as pd
import numpy as np
import collections
import yaml
import copy
import re
from astropy import constants, units as u
from tardis.util import element_symbol2atomic_number

import logging
logger = logging.getLogger(__name__)


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
        A constructor for converting quantity-like YAML strings to
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
        value_str, unit_str = data.split(None, 1)
        value = float(value_str)
        if unit_str.strip() == 'log_lsun':
            value = 10 ** (value + np.log10(constants.L_sun.cgs.value))
            unit_str = 'erg/s'
        return value * u.Unit(unit_str)


YAMLLoader.add_constructor('!quantity', YAMLLoader.construct_quantity)
# This regex matches anything that is a number (scientific notation supported) followed by whitespace
# and any other characters (which should be a unit).
pattern = re.compile(r'^-?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?\s+.+$')
YAMLLoader.add_implicit_resolver('!quantity', pattern)


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
