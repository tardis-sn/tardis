"""
A collection of subpackages and submodules to handle input and output data.
"""

from tardis.io.configuration.config_internal import (
    get_data_dir,
    get_internal_configuration,
)
from tardis.io.model.readers.base import read_density_file

# readin model_data
from tardis.io.model.readers.generic_readers import (
    read_simple_ascii_abundances,
    read_simple_ascii_density,
)
