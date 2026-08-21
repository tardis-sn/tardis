from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt
import pandas as pd
from astropy import units as u

from tardis.io.configuration.config_reader import Configuration

if TYPE_CHECKING:
    from tardis.model.geometry.radial1d import Radial1DGeometry
    from tardis.model.geometry.radial1d_homologous import (
        HomologousRadial1DGeometry,
    )


@dataclass
class CSVYData:
    """
    Data structure for CSVY model data.

    Parameters
    ----------
    model_config : Configuration
        Validated configuration object from the CSVY file.
    velocity : numpy.ndarray
        Velocity array for the model shells.
    density : numpy.ndarray or None
        Density array for the model shells.
    mass_fractions : pandas.DataFrame, optional
        Mass fractions DataFrame with atomic_number as index.
    isotope_mass_fractions : pandas.DataFrame, optional
        Isotope mass fractions DataFrame with MultiIndex of atomic and mass
        number.
    raw_csv_data : pandas.DataFrame or None, optional
        Raw CSV data from the CSVY file.
    """

    model_config: Configuration
    velocity: npt.NDArray[np.float64]
    density: npt.NDArray[np.float64] | None
    mass_fractions: pd.DataFrame = field(default_factory=pd.DataFrame)
    isotope_mass_fractions: pd.DataFrame = field(default_factory=pd.DataFrame)
    raw_csv_data: pd.DataFrame | None = None

    def to_geometry(
        self, time_explosion: u.Quantity | None = None
    ) -> HomologousRadial1DGeometry:
        """
        Construct a homologous radial geometry from this CSVY data.

        Parameters
        ----------
        time_explosion : astropy.units.Quantity, optional
            Time of explosion. If None, attempts to extract from
            ``model_config``.

        Returns
        -------
        tardis.model.geometry.radial1d_homologous.HomologousRadial1DGeometry
            The geometry object constructed from the CSVY data.
        """
        from tardis.model.geometry.radial1d_homologous import (  # noqa: PLC0415
            HomologousRadial1DGeometry,
        )

        if time_explosion is None:
            # Try to extract time_explosion from model_config
            if hasattr(self.model_config, "time_explosion"):
                time_explosion = self.model_config.time_explosion
                if not isinstance(time_explosion, u.Quantity):
                    time_explosion = u.Quantity(time_explosion)

        v_inner = self.velocity[:-1]
        v_outer = self.velocity[1:]
        geometry = HomologousRadial1DGeometry(
            v_inner=v_inner,
            v_outer=v_outer,
            v_inner_boundary=None,
            v_outer_boundary=None,
            time_explosion=time_explosion,
        )
        return geometry

    def to_nonhomologous_geometry(
        self, time_explosion: u.Quantity | None = None
    ) -> Radial1DGeometry:
        """
        Construct a nonhomologous radial geometry from this CSVY data.

        Parameters
        ----------
        time_explosion : astropy.units.Quantity, optional
            Time of explosion. If None, attempts to extract from model_config.

        Returns
        -------
        tardis.model.geometry.radial1d.Radial1DGeometry
            The geometry object constructed from the CSVY data.
        """
        from tardis.model.geometry.radial1d import (  # noqa: PLC0415
            Radial1DGeometry,
        )

        if time_explosion is None:
            # Try to extract time_explosion from model_config
            if hasattr(self.model_config, "time_explosion"):
                time_explosion = self.model_config.time_explosion
                if not isinstance(time_explosion, u.Quantity):
                    time_explosion = u.Quantity(time_explosion)

        v_inner = self.velocity[:-1]
        v_outer = self.velocity[1:]
        geometry = Radial1DGeometry(
            r_inner=(v_inner * time_explosion).cgs,
            r_outer=(v_outer * time_explosion).cgs,
            v_inner=v_inner,
            v_outer=v_outer,
            r_inner_boundary=None,
            r_outer_boundary=None,
            v_inner_boundary=None,
            v_outer_boundary=None,
        )
        return geometry
