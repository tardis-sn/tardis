from dataclasses import dataclass
import numpy as np
import pandas as pd
from astropy import units as u


@dataclass
class SNECIsotopeProfile:
    enclosed_mass: np.ndarray
    radius: np.ndarray
    isotope_mass_fraction: pd.DataFrame


def read_snec_input_profile(file_path: str) -> SNECIsotopeProfile:
    with open(file_path, "r") as file_handle:
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
        data = pd.read_csv(file_handle, sep="\s+", header=None)

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
        ), f"xData dimensions {isotope_mass_fraction.shape} do not match expected ({row_count}, {column_count - 2})"

    return SNECIsotopeProfile(
        enclosed_mass=enclosed_mass,
        radius=radius,
        isotope_mass_fraction=isotope_mass_fraction,
    )
