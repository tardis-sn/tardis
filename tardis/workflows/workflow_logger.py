from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt
import pandas as pd
from IPython.display import display

from tardis.io.configuration.config_reader import Configuration
from tardis.io.logger.logger import logging_state
from tardis.util.environment import Environment

if TYPE_CHECKING:
    from astropy import units as u

logger = logging.getLogger(__name__)


class WorkflowLogger:
    """Configure workflow logging and report radiative-state changes."""

    def __init__(
        self,
        configuration: Configuration,
        log_level: str | None = None,
        specific_log_level: bool | None = None,
    ) -> None:
        """Create a workflow logger and configure TARDIS logging.

        Parameters
        ----------
        configuration : Configuration
            Configuration object containing the logging settings.
        log_level : str or None, optional
            Override for the configured logging level.
        specific_log_level : bool or None, optional
            Override for the configured specific logging setting.
        """
        logging_state(log_level, configuration, specific_log_level)

    def log_plasma_state(
        self,
        t_rad: u.Quantity,
        dilution_factor: npt.ArrayLike,
        t_inner: u.Quantity,
        next_t_rad: u.Quantity,
        next_dilution_factor: npt.ArrayLike,
        next_t_inner: u.Quantity,
        log_sampling: int = 5,
    ) -> None:
        """Log the current and estimated plasma state.

        Parameters
        ----------
        t_rad : astropy.units.Quantity
            Current radiation temperature by shell.
        dilution_factor : numpy.ndarray
            Current dilution factor by shell.
        t_inner : astropy.units.Quantity
            Current inner-boundary temperature.
        next_t_rad : astropy.units.Quantity
            Estimated radiation temperature for the next iteration.
        next_dilution_factor : numpy.ndarray
            Estimated dilution factor for the next iteration.
        next_t_inner : astropy.units.Quantity
            Estimated inner-boundary temperature for the next iteration.
        log_sampling : int, optional
            Log every ``log_sampling``-th shell.
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

        logger.info("\n\tPlasma stratification:")
        if Environment.allows_widget_display():
            if logger.level <= logging.INFO and (
                not logger.filters or logger.filters[0].log_level == logging.INFO
            ):
                display(
                    plasma_state_log.iloc[::log_sampling].style.format("{:.3g}")
                )
        else:
            output = plasma_state_log.iloc[::log_sampling].to_string(
                float_format=lambda value: f"{value:.3g}",
                justify="center",
            )
            output_df = "".join(f"\t{line}\n" for line in output.split("\n"))
            logger.info("\n%s", output_df)

        logger.info(
            "\n\tCurrent t_inner = %s\n"
            "\tExpected t_inner for next iteration = %s\n",
            f"{t_inner:.3f}",
            f"{next_t_inner:.3f}",
        )
