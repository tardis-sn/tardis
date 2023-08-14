"""
A collection of subpackages and submodules to handle input and output data. 
"""

# readin model_data
from tardis.io.model.readers.generic_readers import (
    read_simple_ascii_abundances,
)

from tardis.io.config_internal import get_internal_configuration, get_data_dir
from tardis.io.model.readers.generic_readers import read_density_file, read_simple_ascii_density
