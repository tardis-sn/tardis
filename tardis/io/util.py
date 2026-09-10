# Utility functions for the IO part of TARDIS

import collections.abc as collections_abc
import hashlib
import logging
import shutil
from collections import OrderedDict
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np
import yaml
from astropy import units as u
from astropy.utils.data import download_file

from tardis import __path__ as TARDIS_PATH
from tardis import constants as const

logger = logging.getLogger(__name__)


def get_internal_data_path(fname: str) -> str:
    """
    Get internal data path of TARDIS.

    Parameters
    ----------
    fname : str
        The filename to join with the internal data path.

    Returns
    -------
    str
        Internal data path of TARDIS joined with the filename.
    """
    return str(Path(TARDIS_PATH[0]) / "data" / fname)


def quantity_from_str(text: str) -> u.Quantity:
    """
    Convert a string to `astropy.units.Quantity`.

    Parameters
    ----------
    text : str
        The string to convert to `astropy.units.Quantity`. Expected format
        is "value unit", e.g., "1.0 cm" or "5 log_lsun".

    Returns
    -------
    astropy.units.Quantity
        The converted quantity with appropriate units.

    Notes
    -----
    Special handling for "log_lsun" unit which is converted to solar luminosity
    in CGS units.
    """
    value_str, unit_str = text.split(None, 1)
    value = float(value_str)
    if unit_str.strip() == "log_lsun":
        value = 10 ** (value + np.log10(const.L_sun.cgs.value))
        unit_str = "erg/s"

    unit = u.Unit(unit_str)
    if unit == u.L_sun:
        return value * const.L_sun

    return u.Quantity(value, unit_str)


class MockRegexPattern:
    """
    A mock class to be used in place of a compiled regular expression
    when a type check is needed instead of a regex match.

    Notes
    -----
    This is usually a lot slower than regex matching.

    Parameters
    ----------
    target_type : type
        The target type for conversion testing.
    """

    def __init__(self, target_type: type) -> None:
        """
        Initialize the MockRegexPattern.

        Parameters
        ----------
        target_type : type
            The target type for conversion testing.
        """
        self.type = target_type

    def match(self, text: str) -> bool:
        """
        Test if text can be converted to the target type.

        Parameters
        ----------
        text : str
            A string to be passed to `target_type` for conversion.

        Returns
        -------
        bool
            Returns `True` if `text` can be converted to `target_type`,
            otherwise returns `False`.
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

    def construct_quantity(self, node: yaml.ScalarNode) -> u.Quantity:
        """
        A constructor for converting quantity-like YAML nodes to
        `astropy.units.Quantity` objects.

        Parameters
        ----------
        node : yaml.Node
            The YAML node to be constructed.

        Returns
        -------
        astropy.units.Quantity
            The constructed quantity object.
        """
        data = self.construct_scalar(node)
        return quantity_from_str(data)

    def mapping_constructor(self, node: yaml.MappingNode) -> OrderedDict:
        """
        Construct an OrderedDict from a YAML mapping node.

        Parameters
        ----------
        node : yaml.Node
            The YAML mapping node to construct.

        Returns
        -------
        OrderedDict
            The constructed ordered dictionary.
        """
        return OrderedDict(self.construct_pairs(node))


YAMLLoader.add_constructor("!quantity", YAMLLoader.construct_quantity)
YAMLLoader.add_implicit_resolver(
    "!quantity", MockRegexPattern(quantity_from_str), None
)
YAMLLoader.add_implicit_resolver(
    "tag:yaml.org,2002:float", MockRegexPattern(float), None
)
YAMLLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
    YAMLLoader.mapping_constructor,
)


def yaml_load_file(filename: str, loader: type = yaml.Loader) -> Any:
    """
    Load a YAML file using the specified loader.

    Parameters
    ----------
    filename : str
        Path to the YAML file to load.
    loader : type, optional
        YAML loader class to use, by default yaml.Loader.

    Returns
    -------
    Any
        The loaded YAML content.
    """
    with open(filename) as stream:
        return yaml.load(stream, Loader=loader)


def traverse_configs(base: Any, other: Any, func: Any, *args: Any) -> None:
    """
    Recursively traverse a base dict or list along with another one
    calling `func` for leafs of both objects.

    Parameters
    ----------
    base : Any
        The object on which the traversing is done.
    other : Any
        The object which is traversed along with `base`.
    func : Any
        A function called for each leaf of `base` and the corresponding leaf of `other`.
        Signature: `func(item1, item2, *args)`.
    *args : Any
        Arguments passed into `func`.
    """
    if isinstance(base, collections_abc.Mapping):
        for k in base:
            traverse_configs(base[k], other[k], func, *args)
    elif (
        isinstance(base, collections_abc.Iterable)
        and not isinstance(base, str)
        and not hasattr(base, "shape")
    ):
        for val1, val2 in zip(base, other):
            traverse_configs(val1, val2, func, *args)
    else:
        func(base, other, *args)


def assert_equality(item1: Any, item2: Any) -> None:
    """
    Assert that two items are equal, handling special cases for units and arrays.

    Parameters
    ----------
    item1 : Any
        First item to compare.
    item2 : Any
        Second item to compare.

    Raises
    ------
    AssertionError
        If the items are not equal.
    """
    assert type(item1) is type(item2)
    try:
        if hasattr(item1, "unit"):
            assert item1.unit == item2.unit
        assert np.allclose(item1, item2, atol=0.0)
    except (ValueError, TypeError):
        assert item1 == item2


def check_equality(item1: Any, item2: Any) -> bool:
    """
    Check if two items are equal using traverse_configs and assert_equality.

    Parameters
    ----------
    item1 : Any
        First item to compare.
    item2 : Any
        Second item to compare.

    Returns
    -------
    bool
        True if items are equal, False otherwise.
    """
    try:
        traverse_configs(item1, item2, assert_equality)
    except AssertionError:
        return False
    else:
        return True


@cache
def download_from_url(
    url: str,
    dst: str,
    checksum: str,
    src: tuple[str, ...] | None = None,
    retries: int = 3,
) -> None:
    """
    Download files from a given URL.

    Parameters
    ----------
    url : str
        URL to download from.
    dst : str
        Destination folder for the downloaded file.
    checksum : str
        Expected MD5 checksum of the file.
    src : tuple of str, optional
        List of URLs to use as mirrors.
    retries : int, optional
        Number of retry attempts, by default 3.

    Raises
    ------
    RuntimeError
        If maximum number of retries is reached and checksum still doesn't match.
    """
    cached_file_path = download_file(url, sources=src, pkgname="tardis")

    with open(cached_file_path, "rb") as f:
        new_checksum = hashlib.md5(f.read()).hexdigest()

    if checksum == new_checksum:
        shutil.copy(cached_file_path, dst)

    elif checksum != new_checksum and retries > 0:
        retries -= 1
        logger.warning(
            "Incorrect checksum, retrying... (%d attempts remaining)",
            retries + 1,
        )
        download_from_url(url, dst, checksum, src, retries)

    else:
        logger.error("Maximum number of retries reached. Aborting")
