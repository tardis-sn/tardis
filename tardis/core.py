#functions that are important for the general usage of TARDIS

from tardis.io import config_reader
from tardis import model, simulation


def run_tardis(configuration_dict, atom_data=None):
    """
    This function is one of the core functions to run TARDIS from a given
    config object.

    It will return a model object containing

    Parameters
    ----------

    configuration_dict: ~dict

    atom_data: ~tardis.atomic.AtomData
        Atomic data to use for this TARDIS simulation. If set to None, the
        atomic data will be loaded according to keywords set in the configuration
        [default=None]
    """

    tardis_config = config_reader.Configuration.from_config_dict(
        configuration_dict, atom_data=atom_data)
    radial1d_mdl = model.Radial1DModel(tardis_config)

    simulation.run_radial1d(radial1d_mdl)

    return radial1d_mdl