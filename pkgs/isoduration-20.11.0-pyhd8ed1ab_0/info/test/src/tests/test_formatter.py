import pytest

from isoduration import DurationFormattingException, format_duration
from isoduration.types import DateDuration, Duration, TimeDuration


@pytest.mark.parametrize(
    "duration, duration_str",
    (
        # Zero.
        (Duration(DateDuration(), TimeDuration()), "P0D"),
        # All fields.
        (
            Duration(
                DateDuration(years=3, months=6, days=4),
                TimeDuration(hours=12, minutes=30, seconds=5),
            ),
            "P3Y6M4DT12H30M5S",
        ),
        (
            Duration(
                DateDuration(years=18, months=9, days=4),
                TimeDuration(hours=11, minutes=9, seconds=8),
            ),
            "P18Y9M4DT11H9M8S",
        ),
        # All fields, only date.
        (
            Duration(DateDuration(years=4, months=5, days=18), TimeDuration()),
            "P4Y5M18D",
        ),
        # All fields, only time.
        (
            Duration(DateDuration(), TimeDuration(hours=2, minutes=3, seconds=4)),
            "PT2H3M4S",
        ),
        # Some fields, date only.
        (Duration(DateDuration(years=4), TimeDuration()), "P4Y"),
        (Duration(DateDuration(weeks=2), TimeDuration()), "P2W"),
        (Duration(DateDuration(months=1), TimeDuration()), "P1M"),
        (Duration(DateDuration(days=6), TimeDuration()), "P6D"),
        # Some fields, time only.
        (Duration(DateDuration(), TimeDuration(hours=2)), "PT2H"),
        (Duration(DateDuration(), TimeDuration(hours=36)), "PT36H"),
        (Duration(DateDuration(), TimeDuration(minutes=1)), "PT1M"),
        (Duration(DateDuration(), TimeDuration(seconds=22)), "PT22S"),
        (
            Duration(DateDuration(), TimeDuration(minutes=3, seconds=4)),
            "PT3M4S",
        ),
        (
            Duration(DateDuration(), TimeDuration(hours=6, seconds=59)),
            "PT6H59S",
        ),
        # Some fields, date and time.
        (
            Duration(
                DateDuration(days=1),
                TimeDuration(hours=2, minutes=3, seconds=4),
            ),
            "P1DT2H3M4S",
        ),
        (
            Duration(DateDuration(weeks=3), TimeDuration(hours=2, seconds=59)),
            "P3WT2H59S",
        ),
        (
            Duration(DateDuration(days=1), TimeDuration(hours=2)),
            "P1DT2H",
        ),
        (
            Duration(DateDuration(days=1), TimeDuration(hours=12)),
            "P1DT12H",
        ),
        (
            Duration(DateDuration(days=23), TimeDuration(hours=23)),
            "P23DT23H",
        ),
        (
            Duration(DateDuration(days=1), TimeDuration(hours=2, minutes=3)),
            "P1DT2H3M",
        ),
        # Floating point.
        (Duration(DateDuration(years=0.5), TimeDuration()), "P0.5Y"),
        (
            Duration(DateDuration(), TimeDuration(hours=8.5, seconds=3)),
            "PT8.5H3S",
        ),
        (Duration(DateDuration(), TimeDuration(hours=2.3)), "PT2.3H"),
        (
            Duration(DateDuration(), TimeDuration(seconds=22.22)),
            "PT22.22S",
        ),
        # Scientific notation.
        (Duration(DateDuration(years=1e9), TimeDuration()), "P1e+09Y"),
        (Duration(DateDuration(years=1e-9), TimeDuration()), "P1e-09Y"),
        (Duration(DateDuration(years=-1e9), TimeDuration()), "-P1e+09Y"),
        (Duration(DateDuration(), TimeDuration(seconds=90e9)), "PT9e+10S"),
        (
            Duration(DateDuration(years=1e3), TimeDuration(hours=1e7)),
            "P1000YT1e+07H",
        ),
        (
            Duration(DateDuration(years=1e3), TimeDuration(hours=1e-7)),
            "P1000YT1e-07H",
        ),
        (
            Duration(DateDuration(years=1e3), TimeDuration(hours=-1e7)),
            "P1000YT-1e+07H",
        ),
        # Signs.
        (Duration(DateDuration(years=-2), TimeDuration()), "-P2Y"),
        (Duration(DateDuration(weeks=-2), TimeDuration()), "-P2W"),
        (
            Duration(DateDuration(weeks=-2.2), TimeDuration()),
            "-P2.2W",
        ),
        (
            Duration(
                DateDuration(years=-3, months=-6, days=-4),
                TimeDuration(hours=-12, minutes=-30, seconds=-5),
            ),
            "-P3Y6M4DT12H30M5S",
        ),
        (
            Duration(
                DateDuration(years=-3, months=-6, days=-4),
                TimeDuration(hours=12, minutes=30, seconds=5),
            ),
            "P-3Y-6M-4DT12H30M5S",
        ),
        (
            Duration(
                DateDuration(years=-3, months=-6, days=-4),
                TimeDuration(hours=-12, minutes=30, seconds=-5),
            ),
            "P-3Y-6M-4DT-12H30M-5S",
        ),
        (
            Duration(
                DateDuration(days=-1),
                TimeDuration(hours=-2, minutes=-3, seconds=-4),
            ),
            "-P1DT2H3M4S",
        ),
        (
            Duration(
                DateDuration(years=-3, months=-6, days=-4),
                TimeDuration(),
            ),
            "-P3Y6M4D",
        ),
        (
            Duration(DateDuration(), TimeDuration(hours=-6, minutes=-3)),
            "-PT6H3M",
        ),
        (
            Duration(DateDuration(), TimeDuration(hours=-6, minutes=3)),
            "PT-6H3M",
        ),
        (
            Duration(DateDuration(), TimeDuration(hours=6, minutes=-3)),
            "PT6H-3M",
        ),
    ),
)
def test_format_duration(duration, duration_str):
    assert format_duration(duration) == duration_str


@pytest.mark.parametrize(
    "duration, exception, error_msg",
    (
        (
            Duration(
                DateDuration(years=3, months=6, days=4, weeks=12),
                TimeDuration(),
            ),
            DurationFormattingException,
            r"Weeks are incompatible with other date designators",
        ),
    ),
)
def test_format_duration_errors(duration, exception, error_msg):
    with pytest.raises(exception) as exc:
        format_duration(duration)

    assert exc.match(error_msg)
