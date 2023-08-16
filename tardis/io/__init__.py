"""
A collection of subpackages and submodules to handle input and output data. 
"""

from tardis.io.configuration.config_internal import get_internal_configuration, get_data_dir

# readin model_data

from tardis.io.model.readers.generic_readers import (
    read_simple_ascii_abundances,
)
from tardis.io.model.readers.generic_readers import (
    read_simple_ascii_density,
)
from tardis.io.model.readers.base import read_density_file