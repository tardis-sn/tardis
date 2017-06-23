# tests for the config reader module
import os
import pytest
from numpy.testing import assert_almost_equal
from jsonschema.exceptions import ValidationError

from tardis.io import config_reader
from tardis.io.config_reader import Configuration


def data_path(filename):
    data_dir = os.path.dirname(__file__)
    return os.path.abspath(os.path.join(data_dir, 'data', filename))


def test_convergence_section_parser():
    test_convergence_section = {'type': 'damped',
                                'lock_t_inner_cyles': 1,
                                't_inner_update_exponent': -0.5,
                                'damping_constant': 0.5,
                                'threshold': 0.05,
                                'fraction': 0.8,
                                'hold_iterations': 3,
                                't_rad': {'damping_constant': 1.0}}

    parsed_convergence_section = config_reader.parse_convergence_section(
        test_convergence_section)

    assert_almost_equal(parsed_convergence_section['t_rad']['damping_constant'],
                        1.0)

    assert_almost_equal(parsed_convergence_section['w']['damping_constant'],
                        0.5)


def test_from_config_dict(tardis_config_verysimple):
    conf = Configuration.from_config_dict(tardis_config_verysimple,
                                          validate=True,
                                          config_dirname='test')
    assert conf.config_dirname == 'test'
    assert_almost_equal(conf.spectrum[0].value,
                        tardis_config_verysimple['spectrum'][0].value)
    assert_almost_equal(conf.spectrum[-1].value,
                        tardis_config_verysimple['spectrum'][-1].value)

    tardis_config_verysimple['spectrum'] = 'Invalid'
    with pytest.raises(ValidationError):
        conf = Configuration.from_config_dict(tardis_config_verysimple,
                                              validate=True,
                                              config_dirname='test')

