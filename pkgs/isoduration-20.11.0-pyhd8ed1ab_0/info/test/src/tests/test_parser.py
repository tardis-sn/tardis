from decimal import Decimal

import pytest

from isoduration.parser import exceptions, parse_duration
from isoduration.types import DateDuration, Duration, TimeDuration


@pytest.mark.parametrize(
    "duration, date_duration, time_duration",
    (
        # Zero.
        ("P0D", DateDuration(), TimeDuration()),
        ("PT0S", DateDuration(), TimeDuration()),
        # All fields.
        (
            "P3Y6M4DT12H30M5S",
            DateDuration(years=3, months=6, days=4),
            TimeDuration(hours=12, minutes=30, seconds=5),
        ),
        (
            "P18Y9M4DT11H9M8S",
            DateDuration(years=18, months=9, days=4),
            TimeDuration(hours=11, minutes=9, seconds=8),
        ),
        # All fields, only date.
        ("P4Y5M18D", DateDuration(years=4, months=5, days=18), TimeDuration()),
        # All fields, only time.
        ("PT2H3M4S", DateDuration(), TimeDuration(hours=2, minutes=3, seconds=4)),
        # Some fields, date only.
        ("P4Y", DateDuration(years=4), TimeDuration()),
        ("P2W", DateDuration(weeks=2), TimeDuration()),
        ("P1M", DateDuration(months=1), TimeDuration()),
        ("P6D", DateDuration(days=6), TimeDuration()),
        # Some fields, time only.
        ("PT2H", DateDuration(), TimeDuration(hours=2)),
        ("PT36H", DateDuration(), TimeDuration(hours=36)),
        ("PT1M", DateDuration(), TimeDuration(minutes=1)),
        ("PT22S", DateDuration(), TimeDuration(seconds=22)),
        ("PT3M4S", DateDuration(), TimeDuration(minutes=3, seconds=4)),
        ("PT6H59S", DateDuration(), TimeDuration(hours=6, seconds=59)),
        # Some fields, date and time.
        (
            "P1DT2H3M4S",
            DateDuration(days=1),
            TimeDuration(hours=2, minutes=3, seconds=4),
        ),
        (
            "P3WT2H59S",
            DateDuration(weeks=3),
            TimeDuration(hours=2, seconds=59),
        ),
        ("P1DT2H", DateDuration(days=1), TimeDuration(hours=2)),
        ("P1DT12H", DateDuration(days=1), TimeDuration(hours=12)),
        ("P23DT23H", DateDuration(days=23), TimeDuration(hours=23)),
        (
            "P1DT2H3M",
            DateDuration(days=1),
            TimeDuration(hours=2, minutes=3),
        ),
        # Floating point.
        ("P0.5Y", DateDuration(years=Decimal("0.5")), TimeDuration()),
        ("P0,5Y", DateDuration(years=Decimal("0.5")), TimeDuration()),
        (
            "PT8.5H3S",
            DateDuration(),
            TimeDuration(hours=Decimal("8.5"), seconds=3),
        ),
        (
            "PT8,5H3S",
            DateDuration(),
            TimeDuration(hours=Decimal("8.5"), seconds=3),
        ),
        ("PT2.3H", DateDuration(), TimeDuration(hours=Decimal("2.3"))),
        (
            "PT22.22S",
            DateDuration(),
            TimeDuration(seconds=Decimal("22.22")),
        ),
        # Scientific notation.
        (
            "PT1e3S",
            DateDuration(),
            TimeDuration(seconds=Decimal("1000")),
        ),
        (
            "PT1e-3S",
            DateDuration(),
            TimeDuration(seconds=Decimal("0.001")),
        ),
        (
            "P1E11Y",
            DateDuration(years=Decimal("1e11")),
            TimeDuration(),
        ),
        (
            "P1E+11Y",
            DateDuration(years=Decimal("1e11")),
            TimeDuration(),
        ),
        (
            "P-1E11Y",
            DateDuration(years=Decimal("-1e11")),
            TimeDuration(),
        ),
        (
            "PT1.03e2H",
            DateDuration(),
            TimeDuration(hours=Decimal("103")),
        ),
        (
            "-PT1.03e2H",
            DateDuration(),
            TimeDuration(hours=Decimal("-103")),
        ),
        (
            "P10,8E23W",
            DateDuration(weeks=Decimal("10.8e23")),
            TimeDuration(),
        ),
        # Leading zeroes.
        ("PT000022.22", DateDuration(), TimeDuration()),
        (
            "PT000022.22H",
            DateDuration(),
            TimeDuration(hours=Decimal("22.22")),
        ),
        # Signs.
        ("+P11D", DateDuration(days=11), TimeDuration()),
        ("-P2Y", DateDuration(years=-2), TimeDuration()),
        ("-P2W", DateDuration(weeks=-2), TimeDuration()),
        ("-P2.2W", DateDuration(weeks=Decimal("-2.2")), TimeDuration()),
        (
            "-P3Y6M4DT12H30M5S",
            DateDuration(years=-3, months=-6, days=-4),
            TimeDuration(hours=-12, minutes=-30, seconds=-5),
        ),
        (
            "-P1DT2H3M4S",
            DateDuration(days=-1),
            TimeDuration(hours=-2, minutes=-3, seconds=-4),
        ),
        ("P-20M", DateDuration(months=-20), TimeDuration()),
        ("PT-6H3M", DateDuration(), TimeDuration(hours=-6, minutes=3)),
        ("-PT-6H3M", DateDuration(), TimeDuration(hours=6, minutes=-3)),
        ("-PT-6H+3M", DateDuration(), TimeDuration(hours=6, minutes=-3)),
        # Unconventional numbers, beyond usual boundaries.
        ("P390D", DateDuration(days=390), TimeDuration()),
        ("P20M", DateDuration(months=20), TimeDuration()),
        ("P1000W", DateDuration(weeks=1000), TimeDuration()),
        ("PT72H", DateDuration(), TimeDuration(hours=72)),
        ("PT1000000M", DateDuration(), TimeDuration(minutes=1000000)),
        # Alternative format.
        (
            "P0018-09-04T11:09:08",
            DateDuration(years=18, months=9, days=4),
            TimeDuration(hours=11, minutes=9, seconds=8),
        ),
        (
            "P0003-06-04T12:30:00",
            DateDuration(years=3, months=6, days=4),
            TimeDuration(hours=12, minutes=30, seconds=0),
        ),
        (
            "P00030604T123005",
            DateDuration(years=3, months=6, days=4),
            TimeDuration(hours=12, minutes=30, seconds=5),
        ),
        # Alternative format, with sign.
        (
            "-P0003-06-04T12:30:05",
            DateDuration(years=-3, months=-6, days=-4),
            TimeDuration(hours=-12, minutes=-30, seconds=-5),
        ),
        (
            "-P00030604T123005",
            DateDuration(years=-3, months=-6, days=-4),
            TimeDuration(hours=-12, minutes=-30, seconds=-5),
        ),
    ),
)
def test_parse_duration(duration, date_duration, time_duration):
    assert parse_duration(duration) == Duration(date_duration, time_duration)


@pytest.mark.parametrize(
    "duration, exception, error_msg",
    (
        # EmptyDuration.
        ("", exceptions.EmptyDuration, r"No duration information provided"),
        ("P", exceptions.EmptyDuration, r"No duration information provided"),
        ("1Y2M", exceptions.EmptyDuration, r"No prefix provided"),
        # IncorrectDesignator.
        (
            "P2M1Y",
            exceptions.IncorrectDesignator,
            r"Wrong date designator, or designator in the wrong order: Y",
        ),
        (
            "P1D2H",
            exceptions.IncorrectDesignator,
            r"Wrong date designator, or designator in the wrong order: H",
        ),
        (
            "P4W2D",
            exceptions.IncorrectDesignator,
            r"Wrong date designator, or designator in the wrong order: D",
        ),
        (
            "P4Y2W",
            exceptions.IncorrectDesignator,
            r"Week is incompatible with any other date designator",
        ),
        (
            "P36w",
            exceptions.IncorrectDesignator,
            r"Wrong date designator, or designator in the wrong order: w",
        ),
        (
            "PT3S6M",
            exceptions.IncorrectDesignator,
            r"Wrong time designator, or designator in the wrong order: M",
        ),
        (
            "PT36W",
            exceptions.IncorrectDesignator,
            r"Wrong time designator, or designator in the wrong order: W",
        ),
        (
            "PT11ℵ",
            exceptions.IncorrectDesignator,
            r"Wrong time designator, or designator in the wrong order: ℵ",
        ),
        (
            "PT6H8M43S8S",
            exceptions.IncorrectDesignator,
            r"Wrong time designator, or designator in the wrong order: S",
        ),
        # NoTime.
        ("PT", exceptions.NoTime, r"Wanted time, no time provided"),
        ("P20MT", exceptions.NoTime, r"Wanted time, no time provided"),
        # UnknownToken.
        ("P12'3W", exceptions.UnknownToken, r"Token not recognizable: '"),
        ("P 8W", exceptions.UnknownToken, r"Token not recognizable:  "),
        ("PT(╯°□°)╯︵ ┻━┻", exceptions.UnknownToken, r"Token not recognizable: \("),
        ("P∰", exceptions.UnknownToken, r"Token not recognizable: ∰"),
        ("P💩", exceptions.UnknownToken, r"Token not recognizable: 💩"),
        ("P0018/09/04T11:09:08", exceptions.UnknownToken, r"Token not recognizable: /"),
        # UnparseableValue.
        (
            "P1YM5D",
            exceptions.UnparseableValue,
            r"Value could not be parsed as decimal: ",
        ),
        (
            "PTS",
            exceptions.UnparseableValue,
            r"Value could not be parsed as decimal: ",
        ),
        (
            "P12..80Y",
            exceptions.UnparseableValue,
            r"Value could not be parsed as decimal: 12..80",
        ),
        (
            "PT11.0.42S",
            exceptions.UnparseableValue,
            r"Value could not be parsed as decimal: 11.0.42",
        ),
        (
            "P1234T",
            exceptions.UnparseableValue,
            r"Value could not be parsed as datetime: 1234T",
        ),
        (
            "P00030604T12300500",
            exceptions.UnparseableValue,
            r"Value could not be parsed as datetime: 00030604T12300500",
        ),
        (
            "P00030604T12300",
            exceptions.UnparseableValue,
            r"Value could not be parsed as datetime: 00030604T12300",
        ),
        (
            "P0003060T123005",
            exceptions.UnparseableValue,
            r"Value could not be parsed as datetime: 0003060T123005",
        ),
    ),
)
def test_parse_duration_errors(duration, exception, error_msg):
    with pytest.raises(exception) as exc:
        parse_duration(duration)

    assert exc.match(error_msg)
