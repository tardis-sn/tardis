from decimal import Decimal

import pytest

from isoduration.operations import max_day_in_month, mod2, mod3, quot2, quot3


@pytest.mark.parametrize(
    "op1, op2, result",
    (
        ("-1", "3", "-1"),
        ("0", "3", "0"),
        ("1", "3", "0"),
        ("2", "3", "0"),
        ("3", "3", "1"),
        ("3.123", "3", "1"),
    ),
)
def test_f_quotient_2(op1, op2, result):
    assert quot2(Decimal(op1), Decimal(op2)) == Decimal(result)


@pytest.mark.parametrize(
    "op1, op2, result",
    (
        ("-1", "3", "2"),
        ("0", "3", "0"),
        ("1", "3", "1"),
        ("2", "3", "2"),
        ("3", "3", "0"),
        ("3.123", "3", "0.123"),
    ),
)
def test_modulo_2(op1, op2, result):
    assert mod2(Decimal(op1), Decimal(op2)) == Decimal(result)


@pytest.mark.parametrize(
    "op1, op2, op3, result",
    (
        ("0", "1", "13", "-1"),
        ("1", "1", "13", "0"),
        ("2", "1", "13", "0"),
        ("3", "1", "13", "0"),
        ("4", "1", "13", "0"),
        ("5", "1", "13", "0"),
        ("6", "1", "13", "0"),
        ("7", "1", "13", "0"),
        ("8", "1", "13", "0"),
        ("9", "1", "13", "0"),
        ("10", "1", "13", "0"),
        ("11", "1", "13", "0"),
        ("12", "1", "13", "0"),
        ("13", "1", "13", "1"),
        ("13.123", "1", "13", "1"),
    ),
)
def test_f_quotient_3(op1, op2, op3, result):
    assert quot3(Decimal(op1), Decimal(op2), Decimal(op3)) == Decimal(result)


@pytest.mark.parametrize(
    "op1, op2, op3, result",
    (
        ("0", "1", "13", "12"),
        ("1", "1", "13", "1"),
        ("2", "1", "13", "2"),
        ("3", "1", "13", "3"),
        ("4", "1", "13", "4"),
        ("5", "1", "13", "5"),
        ("6", "1", "13", "6"),
        ("7", "1", "13", "7"),
        ("8", "1", "13", "8"),
        ("9", "1", "13", "9"),
        ("10", "1", "13", "10"),
        ("11", "1", "13", "11"),
        ("12", "1", "13", "12"),
        ("13", "1", "13", "1"),
        ("13.123", "1", "13", "1.123"),
    ),
)
def test_modulo_3(op1, op2, op3, result):
    assert mod3(Decimal(op1), Decimal(op2), Decimal(op3)) == Decimal(result)


@pytest.mark.parametrize(
    "year, month, max_day",
    (
        (1987, 1, 31),
        (1987, 2, 28),
        (1987, 3, 31),
        (1987, 4, 30),
        (1987, 5, 31),
        (1987, 6, 30),
        (1987, 7, 31),
        (1987, 8, 31),
        (1987, 9, 30),
        (1987, 10, 31),
        (1987, 11, 30),
        (1987, 12, 31),
        (1987, 12, 31),
        (2000, 2, 29),
        (1900, 2, 28),
        (2004, 2, 29),
    ),
)
def test_max_day_in_month(year, month, max_day):
    assert max_day_in_month(Decimal(year), Decimal(month)) == Decimal(max_day)
