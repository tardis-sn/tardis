from dataclasses import dataclass, field

import pandas as pd
from astropy import units as u

from tardis.model.geometry.radial1d import (
    Radial1DGeometry,
)
from tardis.model.geometry.radial1d_homologous import (
    HomologousRadial1DGeometry,
)


@dataclass
class ArtisData:
    """
    Data structure for ARTIS model data.

    Parameters
    ----------
    time_of_model : astropy.units.Quantity
        Time at which the ARTIS model is valid.
    velocity : astropy.units.Quantity
        Velocity boundaries for the model shells.
    mean_density : astropy.units.Quantity
        Mean density for each model shell.
    mass_fractions : pandas.DataFrame, optional
        Elemental mass fractions indexed by atomic number.
    isotope_mass_fractions : pandas.DataFrame, optional
        Isotopic mass fractions indexed by atomic and mass number.
    """

    time_of_model: u.Quantity
    velocity: u.Quantity
    mean_density: u.Quantity
    mass_fractions: pd.DataFrame = field(default_factory=pd.DataFrame)
    isotope_mass_fractions: pd.DataFrame = field(default_factory=pd.DataFrame)

    def to_geometry(self) -> HomologousRadial1DGeometry:
        """
        Construct a homologous radial geometry from this ARTIS data.

        ``time_of_model`` is used as ``time_explosion``.

        Returns
        -------
        tardis.model.geometry.radial1d_homologous.HomologousRadial1DGeometry
            The geometry object constructed from the ARTIS data.
        """
        v_inner = self.velocity[:-1]
        v_outer = self.velocity[1:]
        geometry = HomologousRadial1DGeometry(
            v_inner=v_inner,
            v_outer=v_outer,
            v_inner_boundary=None,
            v_outer_boundary=None,
            time_explosion=self.time_of_model,
        )
        return geometry

    def to_nonhomologous_geometry(self) -> Radial1DGeometry:
        """
        Construct a nonhomologous radial geometry from this ARTIS data.

        ``time_of_model`` is used to construct the radius boundaries.

        Returns
        -------
        tardis.model.geometry.radial1d.Radial1DGeometry
            The geometry object constructed from the ARTIS data.
        """
        v_inner = self.velocity[:-1]
        v_outer = self.velocity[1:]
        geometry = Radial1DGeometry(
            r_inner=(v_inner * self.time_of_model).cgs,
            r_outer=(v_outer * self.time_of_model).cgs,
            v_inner=v_inner,
            v_outer=v_outer,
            r_inner_boundary=None,
            r_outer_boundary=None,
            v_inner_boundary=None,
            v_outer_boundary=None,
        )
        return geometry
