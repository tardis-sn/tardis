from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from astropy import units as u

from tardis.model.geometry.radial1d import HomologousRadial1DGeometry


@dataclass
class ArtisData:
    time_of_model: u.Quantity
    velocity: np.ndarray
    mean_density: u.Quantity
    mass_fractions: pd.DataFrame = field(default_factory=pd.DataFrame)

    def to_geometry(self) -> HomologousRadial1DGeometry:
        """
        Construct a HomologousRadial1DGeometry object from this ArtisData.

        The time_of_model is used as the time_explosion.

        Returns
        -------
        tardis.model.geometry.radial1d.HomologousRadial1DGeometry
            The geometry object constructed from the ARTIS data.
        """
        geometry = HomologousRadial1DGeometry(
            v_inner=self.velocity[:-1],
            v_outer=self.velocity[1:],
            v_inner_boundary=None,
            v_outer_boundary=None,
            time_explosion=self.time_of_model,
        )
        return geometry
