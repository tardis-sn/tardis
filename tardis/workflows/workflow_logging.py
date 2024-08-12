import logging

import numpy as np
import pandas as pd
from IPython.display import display

from tardis.io.logger.logger import logging_state
from tardis.util.base import is_notebook

logger = logging.getLogger(__name__)


class WorkflowLogging:
    def __init__(
        self,
        configuration,
        log_level=None,
        specific_log_level=None,
    ):
        logging_state(log_level, configuration, specific_log_level)

    def log_iteration_results(
        self, emitted_luminosity, absorbed_luminosity, luminosity_requested
    ):
        """Print current iteration information to log at INFO level

        Parameters
        ----------
        emitted_luminosity : Quantity
            Current iteration emitted luminosity
        absorbed_luminosity : Quantity
            Current iteration absorbed luminosity
        luminosity_requested : Quantity
            The requested luminosity for the simulation
        """
        logger.info(
            f"\n\tLuminosity emitted   = {emitted_luminosity:.3e}\n"
            f"\tLuminosity absorbed  = {absorbed_luminosity:.3e}\n"
            f"\tLuminosity requested = {luminosity_requested:.3e}\n"
        )

    def log_plasma_state(
        self,
        t_rad,
        dilution_factor,
        t_inner,
        next_t_rad,
        next_dilution_factor,
        next_t_inner,
        log_sampling=5,
    ):
        """
        Logging the change of the plasma state

        Parameters
        ----------
        t_rad : astropy.units.Quanity
            current t_rad
        dilution_factor : np.ndarray
            current dilution_factor
        next_t_rad : astropy.units.Quanity
            next t_rad
        next_dilution_factor : np.ndarray
            next dilution_factor
        log_sampling : int
            the n-th shells to be plotted

        Returns
        -------
        """
        plasma_state_log = pd.DataFrame(
            index=np.arange(len(t_rad)),
            columns=["t_rad", "next_t_rad", "w", "next_w"],
        )
        plasma_state_log["t_rad"] = t_rad
        plasma_state_log["next_t_rad"] = next_t_rad
        plasma_state_log["w"] = dilution_factor
        plasma_state_log["next_w"] = next_dilution_factor
        plasma_state_log.columns.name = "Shell No."

        if is_notebook():
            logger.info("\n\tPlasma stratification:")

            # Displaying the DataFrame only when the logging level is NOTSET, DEBUG or INFO
            if logger.level <= logging.INFO:
                if not logger.filters:
                    display(
                        plasma_state_log.iloc[::log_sampling].style.format(
                            "{:.3g}"
                        )
                    )
                elif logger.filters[0].log_level == 20:
                    display(
                        plasma_state_log.iloc[::log_sampling].style.format(
                            "{:.3g}"
                        )
                    )
        else:
            output_df = ""
            plasma_output = plasma_state_log.iloc[::log_sampling].to_string(
                float_format=lambda x: f"{x:.3g}",
                justify="center",
            )
            for value in plasma_output.split("\n"):
                output_df = output_df + f"\t{value}\n"
            logger.info("\n\tPlasma stratification:")
            logger.info(f"\n{output_df}")

        logger.info(
            f"\n\tCurrent t_inner = {t_inner:.3f}\n\tExpected t_inner for next iteration = {next_t_inner:.3f}\n"
        )
