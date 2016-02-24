#functions that are important for the general usage of TARDIS

def run_tardis(config, atom_data=None):
    """
    This function is one of the core functions to run TARDIS from a given
    config object.

    It will return a model object containing

    Parameters
    ----------

    config: ~str or ~dict
        filename of configuration yaml file or dictionary

    atom_data: ~str or ~tardis.atomic.AtomData
        if atom_data is a string it is interpreted as a path to a file storing
        the atomic data. Atomic data to use for this TARDIS simulation. If set to None, the
        atomic data will be loaded according to keywords set in the configuration
        [default=None]
    """
    import yaml

    from tardis.io import config_reader
    from tardis import model, simulation, atomic

    try:
        config_dict = yaml.load(open(config))
    except TypeError:
        config_dict = config

    if atom_data is not None:
        try:
            atom_data = atomic.AtomData.from_hdf5(atom_data)
        except TypeError:
            atom_data = atom_data

    tardis_config = config_reader.Configuration.from_config_dict(
        config_dict, atom_data=atom_data)
    radial1d_mdl = model.Radial1DModel(tardis_config)

    simulation.run_radial1d(radial1d_mdl)

    return radial1d_mdl