#testing suite for config reader


from tardis import config_reader
from astropy import units as u
import os

default_config = os.path.join(os.path.dirname(config_reader.__file__), 'data', 'tardis_example_config.yml')

tardis_config = config_reader.TARDISConfiguration.from_yaml(default_config)


def test_time_explosion():
    assert tardis_config.t_explosion == (13 * u.day).to('s').value