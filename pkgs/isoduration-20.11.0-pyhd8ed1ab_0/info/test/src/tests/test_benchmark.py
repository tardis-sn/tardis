from isodate import parse_duration
from isoduration import parse_duration as isodate_parse_duration


def test_isoduration(benchmark):
    benchmark(parse_duration, "P18Y9M4DT11H9M8S")


def test_isodate(benchmark):
    benchmark(isodate_parse_duration, "P18Y9M4DT11H9M8S")
