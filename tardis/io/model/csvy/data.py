from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from astropy import units as u

from tardis.io.configuration.config_reader import Configuration


@dataclass
class CSVYData:
    """
    Data structure for CSVY model data.

    Parameters
    ----------
    model_config : Configuration
        Validated configuration object from the CSVY file.
    velocity : np.ndarray
        Velocity array for the model shells.
    density : np.ndarray or None
        Density array for the model shells.
    mass_fractions : pd.DataFrame, optional
        Mass fractions DataFrame with atomic_number as index.
    isotope_mass_fractions : pd.DataFrame, optional
        Isotope mass fractions DataFrame with MultiIndex of (atomic_number, mass_number).
    raw_csv_data : pd.DataFrame or None, optional
        Raw CSV data from the CSVY file.
    """

    model_config: Configuration
    velocity: np.ndarray
    density: np.ndarray | None
    mass_fractions: pd.DataFrame = field(default_factory=pd.DataFrame)
    isotope_mass_fractions: pd.DataFrame = field(default_factory=pd.DataFrame)
    raw_csv_data: pd.DataFrame | None = None

    def to_geometry(self, time_explosion: u.Quantity | None = None):
        """
        Construct a HomologousRadial1DGeometry object from this CSVYData.

        Parameters
        ----------
        time_explosion : astropy.units.Quantity, optional
            Time of explosion. If None, attempts to extract from model_config.

        Returns
        -------
        HomologousRadial1DGeometry
            The geometry object constructed from the CSVY data.
        """
        from tardis.model.geometry.radial1d import HomologousRadial1DGeometry

        if time_explosion is None:
            # Try to extract time_explosion from model_config
            if hasattr(self.model_config, "time_explosion"):
                time_explosion = self.model_config.time_explosion
                if not isinstance(time_explosion, u.Quantity):
                    time_explosion = u.Quantity(time_explosion)

        geometry = HomologousRadial1DGeometry(
            v_inner=self.velocity[:-1],
            v_outer=self.velocity[1:],
            v_inner_boundary=None,
            v_outer_boundary=None,
            time_explosion=time_explosion,
        )
        return geometry
