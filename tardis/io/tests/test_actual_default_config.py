# This is too test the actual default config

import tardis
from tardis.io.default_config_parser import DefaultParser
import os
import yaml
default_config_fname = os.path.join(tardis.__path__[0], 'data', 'tardis_default_config_definition.yml')


def test_read():
    default_config_dict = yaml.load(file(default_config_fname))
    DefaultParser(default_config_dict)