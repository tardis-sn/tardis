from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from astropy import units as u

from tardis.model.geometry.radial1d import HomologousRadial1DGeometry


@dataclass
class CSVYData:
    """
    Data structure for CSVY model data.

    Parameters
    ----------
    yaml_dict : dict
        YAML metadata from the CSVY file.
    velocity : np.ndarray
        Velocity array for the model shells.
    density : np.ndarray
        Density array for the model shells.
    mass_fractions : pd.DataFrame, optional
        Mass fractions DataFrame with atomic_number as index.
    isotope_mass_fractions : pd.DataFrame, optional
        Isotope mass fractions DataFrame with MultiIndex of (atomic_number, mass_number).
    """

    yaml_dict: dict
    velocity: np.ndarray
    density: np.ndarray
    mass_fractions: pd.DataFrame = field(default_factory=pd.DataFrame)
    isotope_mass_fractions: pd.DataFrame = field(default_factory=pd.DataFrame)

    def to_geometry(
        self, time_explosion: u.Quantity | None = None
    ) -> HomologousRadial1DGeometry:
        """
        Construct a HomologousRadial1DGeometry object from this CSVYData.

        Parameters
        ----------
        time_explosion : astropy.units.Quantity, optional
            Time of explosion. If None, attempts to extract from yaml_dict.

        Returns
        -------
        tardis.model.geometry.radial1d.HomologousRadial1DGeometry
            The geometry object constructed from the CSVY data.
        """
        if time_explosion is None:
            # Try to extract time_explosion from yaml_dict
            time_explosion = self.yaml_dict.get("time_explosion")
            if time_explosion is not None and not isinstance(
                time_explosion, u.Quantity
            ):
                time_explosion = u.Quantity(time_explosion)

        geometry = HomologousRadial1DGeometry(
            v_inner=self.velocity[:-1],
            v_outer=self.velocity[1:],
            v_inner_boundary=None,
            v_outer_boundary=None,
            time_explosion=time_explosion,
        )
        return geometry
