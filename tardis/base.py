"""
Functions that are important for the general usage of TARDIS.
"""

import logging


logger = logging.getLogger(__name__)


def run_tardis(
    config,
    atom_data=None,
    packet_source=None,
    simulation_callbacks=[],
    virtual_packet_logging=False,
    show_convergence_plots=False,
    log_level=None,
    specific_log_level=None,
    show_progress_bars=True,
    **kwargs,
):
    """
    Run TARDIS from a given config object.

    It will return a model object containing the TARDIS Simulation.

    Parameters
    ----------
    config : str or dict or tardis.io.config_reader.Configuration
        filename of configuration yaml file or dictionary or TARDIS Configuration object
    atom_data : str or tardis.atomic.AtomData, optional
        If atom_data is a string it is interpreted as a path to a file storing
        the atomic data. Atomic data to use for this TARDIS simulation. If set to None (i.e. default),
        the atomic data will be loaded according to keywords set in the configuration
    packet_source : class, optional
        A custom packet source class or a child class of `tardis.montecarlo.packet_source`
        used to override the TARDIS `BasePacketSource` class.
    simulation_callbacks : list of lists, default: `[]`, optional
        Set of callbacks to call at the end of every iteration of the Simulation.
        The format of the lists should look like:
        [[callback1, callback_arg1], [callback2, callback_arg2], ...],
        where the callback function signature should look like:
        callback_function(simulation, extra_arg1, ...)
    virtual_packet_logging : bool, default: False, optional
        Option to enable virtual packet logging.
    log_level : {'NOTSET', 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'}, default: None, optional
        Set the level of the TARDIS logger (follows native python logging framework log levels).
        Use this parameter to override the `log_level` specified in the configuration file.
        The default value `None` means that the `log_level` specified in the configuration file will be used.
    specific_log_level : bool, default: None, optional
        Allows to set specific logging levels, overriding the value in the configuration file.
        If True, only show the log messages from a particular log level, set by `log_level`.
        If False, the logger shows log messages belonging to the level set and all levels above it in severity.
        The default value None means that the `specific_log_level` specified in the configuration file will be used.
    show_convergence_plots : bool, default: False, optional
        Option to enable tardis convergence plots.
    show_progress_bars : bool, default: True, optional
        Option to enable the progress bar.
    **kwargs : dict, optional
        Optional keyword arguments including those
        supported by :obj:`tardis.visualization.tools.convergence_plot.ConvergencePlots`.


    Returns
    -------
    tardis.simulation.Simulation

    Notes
    -----
    Please see the `logging tutorial <https://tardis-sn.github.io/tardis/io/optional/logging_configuration.html>`_ to know more about `log_level` and `specific` options.
    """
    from tardis.io.logger.logger import logging_state
    from tardis.io.config_reader import Configuration
    from tardis.io.atom_data.base import AtomData
    from tardis.simulation import Simulation

    if isinstance(config, Configuration):
        tardis_config = config
    else:
        try:
            tardis_config = Configuration.from_yaml(config)
        except TypeError:
            logger.debug(
                "TARDIS Config not available via YAML. Reading through TARDIS Config Dictionary"
            )
            tardis_config = Configuration.from_config_dict(config)

    if not isinstance(show_convergence_plots, bool):
        raise TypeError("Expected bool in show_convergence_plots argument")

    logging_state(log_level, tardis_config, specific_log_level)

    if atom_data is not None:
        try:
            atom_data = AtomData.from_hdf(atom_data)
        except TypeError:
            logger.debug(
                "Atom Data Cannot be Read from HDF. Setting to Default Atom Data"
            )
            atom_data = atom_data

    simulation = Simulation.from_config(
        tardis_config,
        packet_source=packet_source,
        atom_data=atom_data,
        virtual_packet_logging=virtual_packet_logging,
        show_convergence_plots=show_convergence_plots,
        show_progress_bars=show_progress_bars,
        **kwargs,
    )
    for cb in simulation_callbacks:
        simulation.add_callback(*cb)

    simulation.run()

    return simulation
