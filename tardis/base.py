"""
Functions that are important for the general usage of TARDIS.
"""

import logging

logger = logging.getLogger(__name__)


def run_tardis(
    config,
    simulation_callbacks=None,
    virtual_packet_logging=False,
    show_convergence_plots=False,
    log_level=None,
    specific_log_level=None,
    show_progress_bars=True,
    display_logging_widget=True,
    **kwargs,
):
    """
    Run TARDIS from a given config object.

    It will return a model object containing the TARDIS Simulation.

    Parameters
    ----------
    config : str or dict or tardis.io.config_reader.Configuration
        filename of configuration yaml file or dictionary or TARDIS Configuration object
    packet_source : class, optional
        A custom packet source class or a child class of `tardis.transport.montecarlo.packet_source`
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
    display_widget : bool, default: True, optional
        Option to display the logging widget in Jupyter/VSCode environments.
        If False, logs will be printed normally instead.
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
    from tardis.io.configuration.config_reader import Configuration
    from tardis.io.logger.logger import logging_state
    from tardis.workflows.standard_tardis_workflow import StandardTARDISWorkflow
    from tardis.io.atom_data.util import resolve_atom_data_fname
    from tardis.io.atom_data import download_atom_data
    from tardis.io.atom_data import AtomData

    if simulation_callbacks is None:
        simulation_callbacks = []
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

    logger_widget, tardislogger = logging_state(log_level, tardis_config, specific_log_level, display_logging_widget)

    convergence_plots_config_options = [
        "plasma_plot_config",
        "t_inner_luminosities_config",
        "plasma_cmap",
        "t_inner_luminosities_colors",
        "export_convergence_plots",
    ]
    convergence_plots_kwargs = {}
    for item in set(convergence_plots_config_options).intersection(kwargs.keys()):
        convergence_plots_kwargs[item] = kwargs[item]

    workflow = StandardTARDISWorkflow(
        configuration=tardis_config,
        enable_virtual_packet_logging=virtual_packet_logging,
        log_level=log_level,
        specific_log_level=specific_log_level,
        show_progress_bars=show_progress_bars,
        show_convergence_plots=show_convergence_plots,
        convergence_plots_kwargs=convergence_plots_kwargs,
    )

    workflow.run()
    if logger_widget:
        tardislogger.finalize_widget_logging()
        tardislogger.remove_widget_handler()
        tardislogger.setup_stream_handler()
    return workflow
