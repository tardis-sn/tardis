# tests for the config reader module

from tardis import config_reader
from astropy import units as u
import pytest



def test_quantity_parser():
    q1 = config_reader.parse2quantity('5 km/s')
    assert q1.value == 5.
    assert q1.unit == u.Unit('km/s')


