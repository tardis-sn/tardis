# functions that are important for the general usage of TARDIS


def run_tardis(
    config,
    atom_data=None,
    packet_source=None,
    simulation_callbacks=[],
    verbosity=None,
    save=None,
):
    """
    This function is one of the core functions to run TARDIS from a given
    config object.

    It will return a model object containing

    Parameters
    ----------
    config : str or dict
        filename of configuration yaml file or dictionary
    atom_data : str or tardis.atomic.AtomData
        if atom_data is a string it is interpreted as a path to a file storing
        the atomic data. Atomic data to use for this TARDIS simulation. If set to None, the
        atomic data will be loaded according to keywords set in the configuration
        [default=None]
    """
    from tardis.io.config_reader import Configuration
    from tardis.io.atom_data.base import AtomData
    from tardis.simulation import Simulation

    import warnings
    import tardis.util.custom_logger as custom_logger_settings

    if verbosity is not None or save is not None:
        if verbosity in ["TARDIS INFO", "ERROR", "CRITICAL"]:
            warnings.filterwarnings("ignore", category=RuntimeWarning)
        custom_logger_settings.init()
        custom_logger_settings.logger.remove()
        if save is not None:
            custom_logger_settings.save = save
        custom_logger_settings.level = verbosity
        custom_logger_settings.reset_logger()

    if atom_data is not None:
        try:
            atom_data = AtomData.from_hdf(atom_data)
        except TypeError:
            atom_data = atom_data

    try:
        tardis_config = Configuration.from_yaml(config)
    except TypeError:
        tardis_config = Configuration.from_config_dict(config)

    simulation = Simulation.from_config(
        tardis_config, packet_source=packet_source, atom_data=atom_data
    )
    for cb in simulation_callbacks:
        simulation.add_callback(*cb)

    simulation.run()

    return simulation
