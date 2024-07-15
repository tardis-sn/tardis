#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `rfc3339_validator` package."""
import re
import pytest
from rfc3339_validator import validate_rfc3339
import strict_rfc3339
from hypothesis import given, settings, example
import hypothesis.strategies as st
import six

# It is supposed to be used to generate both valid and invalid dates
RFC3339_REGEX = r"""
    ^
    (\d\d\d\d) # Year
    -
    (\d\d)     # Month
    -
    (\d\d)     # Day
    [T ]
    (\d\d)     # Hours
    :
    (?:\d\d)     # Minutes
    :
    (?:\d\d)     # Seconds
    (?:\.\d+)?   # Secfrac
    (  [Zz]
      |(?:[+\-])(\d\d):(?:\d\d)
    )
    $
"""
if six.PY3:
    RFC3339_REGEX_FLAG = re.X | re.A
else:
    RFC3339_REGEX_FLAG = re.X
RFC3339_REGEX_ASCII = re.compile(RFC3339_REGEX, RFC3339_REGEX_FLAG)
RFC3339_REGEX_UNICODE = re.compile(RFC3339_REGEX, re.X)


@pytest.mark.skipif(six.PY2, reason="Requires python3 or higher, because strftime on python 2 only supports dates "
                                    "newer than 1900")
@given(datetime_str=st.datetimes().filter(lambda d: d.year > 1000).map(lambda d: d.strftime("%Y-%m-%dT%H:%M:%SZ")))
def test_valid_dates(datetime_str):
    assert validate_rfc3339(datetime_str)


@settings(max_examples=1500)
@given(datetime_str=st.from_regex(RFC3339_REGEX_ASCII, fullmatch=True))
@example(datetime_str='')
def test_against_legacy(datetime_str):
    legacy_result = strict_rfc3339.validate_rfc3339(datetime_str)
    new_result = validate_rfc3339(datetime_str)
    assert legacy_result == new_result


@settings(max_examples=1)
@given(datetime_str=st.from_regex(RFC3339_REGEX_UNICODE, fullmatch=True))
@example(datetime_str='0001-01-01T00:00:0٠Z')
def test_with_unicode(datetime_str):
    assert not validate_rfc3339(datetime_str)
