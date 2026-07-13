from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from astropy import units as u

from tardis.model.geometry.radial1d_nonhomologous import (
    NonhomologousRadial1DGeometry,
)


@dataclass
class ArtisData:
    time_of_model: u.Quantity
    velocity: np.ndarray
    mean_density: u.Quantity
    mass_fractions: pd.DataFrame = field(default_factory=pd.DataFrame)

    def to_geometry(self) -> NonhomologousRadial1DGeometry:
        """
        Construct a NonhomologousRadial1DGeometry object from this ArtisData.

        The time_of_model is used as the time_explosion.

        Returns
        -------
        tardis.model.geometry.radial1d_nonhomologous.NonhomologousRadial1DGeometry
            The geometry object constructed from the ARTIS data.
        """
        v_inner = self.velocity[:-1]
        v_outer = self.velocity[1:]
        geometry = NonhomologousRadial1DGeometry(
            r_inner=(v_inner * self.time_of_model).cgs,
            r_outer=(v_outer * self.time_of_model).cgs,
            v_inner=self.velocity[:-1],
            v_outer=self.velocity[1:],
            r_inner_boundary=None,
            r_outer_boundary=None,
            v_inner_boundary=None,
            v_outer_boundary=None,
        )
        return geometry
