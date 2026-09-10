from dataclasses import dataclass

import numpy as np
import pandas as pd
import xarray as xr
from astropy import units as u


@dataclass
class SNECIsotopeProfile:
    """
    Data class representing an isotope profile from SNEC output.

    Attributes
    ----------
    enclosed_mass : np.ndarray
        Array of enclosed mass values (with astropy units).
    radius : np.ndarray
        Array of radius values (with astropy units).
    isotope_mass_fraction : pd.DataFrame
        DataFrame of isotope mass fractions, indexed by cell and columns as isotopes.
    """

    enclosed_mass: np.ndarray
    radius: np.ndarray
    isotope_mass_fraction: pd.DataFrame

    def to_xr_array(self) -> xr.DataArray:
        """
        Convert the isotope mass fractions into a 2D DataArray with dimensions:
            - cell_id
            - isotope
        and attach radius and enclosed_mass as coordinates on cell_id.

        Returns
        -------
        xr.DataArray
            2D DataArray of isotope mass fractions with coordinates for radius and enclosed mass.
        """
        # Ensure DataFrame index and columns are named
        df = self.isotope_mass_fraction.copy()
        df.index.name = 'cell_id'
        df.columns.name = 'isotope'

        # Build the DataArray for isotope mass fractions
        da = xr.DataArray(
            data=df.values,
            dims=('cell_id', 'isotope'),
            coords={
                'cell_id': df.index,
                'isotope': df.columns,
            },
            name='isotope_mass_fraction'
        )

        # Attach radius and enclosed_mass as coordinates on cell_id
        da = da.assign_coords(
            radius=('cell_id', self.radius.value),
            enclosed_mass=('cell_id', self.enclosed_mass.value)
        )

        return da


def read_snec_isotope_profile(file_path: str) -> SNECIsotopeProfile:
    """
    Read a SNEC isotope profile file and return a SNECIsotopeProfile object.

    Parameters
    ----------
    file_path : str
        Path to the SNEC isotope profile file.

    Returns
    -------
    SNECIsotopeProfile
        Parsed isotope profile data.
    """
    with open(file_path) as file_handle:
        # Read the first line to get row and column counts
        first_line = file_handle.readline().strip()
        row_count, column_count = map(int, first_line.split())

        # Read the second and third lines as integer arrays
        second_line = file_handle.readline().strip().replace("d", "e")
        mass_number = np.array(list(map(float, second_line.split()))).astype(int)
        third_line = file_handle.readline().strip().replace("d", "e")
        neutron_number = np.array(list(map(float, third_line.split()))).astype(int)
        element_number = mass_number - neutron_number

        # Read the remaining lines into a DataFrame
        data = pd.read_csv(file_handle, sep=r"\s+", header=None)

        enclosed_mass = data.iloc[:, 0].values * u.g
        radius = data.iloc[:, 1].values * u.cm

        isotope_mass_fraction = data.iloc[:, 2:]
        isotope_mass_fraction.columns = pd.MultiIndex.from_arrays(
            [element_number, mass_number], names=["element_number", "mass_number"]
        )

        # Assert the dimensions match the provided row and column counts
        assert isotope_mass_fraction.shape == (
            row_count,
            column_count,
        ), f"Data dimensions {isotope_mass_fraction.shape} do not match expected ({row_count}, {column_count - 2})"

    return SNECIsotopeProfile(
        enclosed_mass=enclosed_mass,
        radius=radius,
        isotope_mass_fraction=isotope_mass_fraction,
    )
