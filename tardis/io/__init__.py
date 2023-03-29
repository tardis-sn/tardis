"""
A collection of subpackages and submodules to handle input and output data. 
"""

# readin model_data
from tardis.io.model_reader import (
    read_simple_ascii_density,
    read_simple_ascii_abundances,
    read_density_file,
)

from tardis.io.config_internal import get_internal_configuration, get_data_dir
