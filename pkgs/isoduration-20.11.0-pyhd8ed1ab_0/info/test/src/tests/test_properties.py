from datetime import timedelta

from hypothesis import given
from hypothesis.strategies import SearchStrategy, builds, datetimes, decimals

from isoduration.formatter import format_duration
from isoduration.parser import parse_duration
from isoduration.types import DateDuration, Duration, TimeDuration


item_st = decimals(
    min_value=-1_000_000_000_000, max_value=+1_000_000_000_000, places=10
)

date_duration_st: SearchStrategy[DateDuration] = builds(
    DateDuration, years=item_st, months=item_st, days=item_st
)
time_duration_st: SearchStrategy[TimeDuration] = builds(
    TimeDuration, hours=item_st, minutes=item_st, seconds=item_st
)


@given(date_duration=date_duration_st, time_duration=time_duration_st)
def test_parse_inverse_of_format(date_duration, time_duration):
    duration = Duration(date_duration, time_duration)
    assert parse_duration(format_duration(duration)) == duration


@given(date_duration=date_duration_st, time_duration=time_duration_st)
def test_duration_double_negation(date_duration, time_duration):
    duration = Duration(date_duration, time_duration)
    neg_duration = -duration

    assert -neg_duration == duration


@given(datetimes())
def test_duration_addition_identity_element(base_datetime):
    identity = Duration(DateDuration(), TimeDuration())
    assert (base_datetime + identity) - base_datetime < timedelta(seconds=1)
