"""Writers for CSVY model files."""

from pathlib import Path

import pandas as pd
import yaml
from astropy import units as u

from tardis.configuration.sorting_globals import SORTING_ALGORITHM

YAML_DELIMITER = "---"


def write_csvy(
    path: str | Path,
    density_t_0: u.Quantity,
    isotope_t_0: u.Quantity,
    velocity: u.Quantity,
    density: u.Quantity,
    abundance: pd.DataFrame,
    elements: list[str],
    overwrite: bool = False,
) -> Path:
    """Write a CSVY model file.

    Parameters
    ----------
    path : str or pathlib.Path
        Output path. The ``.csvy`` suffix is added or replaced.
    density_t_0 : astropy.units.Quantity
        Initial time for the density model.
    isotope_t_0 : astropy.units.Quantity
        Initial time for isotope decay.
    velocity : astropy.units.Quantity
        Shell boundary velocities.
    density : astropy.units.Quantity
        Shell densities.
    abundance : pandas.DataFrame
        Abundance values indexed by shell and species.
    elements : list of str
        Species names corresponding to the abundance columns.
    overwrite : bool, optional
        Whether to overwrite an existing file. Defaults to False.

    Returns
    -------
    pathlib.Path
        The path of the written CSVY file.

    Raises
    ------
    FileExistsError
        If the output exists and ``overwrite`` is False.
    """
    path = Path(path).with_suffix(".csvy")
    if path.exists() and not overwrite:
        raise FileExistsError(
            "The file already exists. Click the 'overwrite' checkbox to overwrite it."
        )

    metadata = {
        "model_density_time_0": density_t_0,
        "model_isotope_time_0": isotope_t_0,
        "tardis_model_config_version": "v1.0",
        "datatype": {
            "fields": [
                {"name": "velocity", "unit": "km/s"},
                {"name": "density", "unit": "g/cm^3"},
                *[
                    {"name": element, "desc": f"fractional {element} abundance"}
                    for element in elements
                ],
            ]
        },
        "v_inner_boundary": velocity[0],
        "v_outer_boundary": velocity[-1],
    }

    with path.open("w") as file:
        yaml_output = yaml.dump(metadata, sort_keys=False)
        yaml_lines = yaml_output.split("\n")
        yaml_lines[0] = yaml_lines[-1] = YAML_DELIMITER
        file.write("\n".join(yaml_lines) + "\n")

        data = abundance.T.copy()
        data.columns = elements
        data.loc[-1] = [0] * len(elements)
        data.index += 1
        data = data.sort_index(kind=SORTING_ALGORITHM)

        formatted_velocity = pd.Series(
            [f"{value:.3e}" for value in velocity.value], dtype=object
        )
        formatted_velocity.iloc[0] = velocity.value[0]
        formatted_velocity.iloc[-1] = velocity.value[-1]
        data.insert(0, "velocity", formatted_velocity)
        data.insert(1, "density", pd.Series(density.value))
        data.to_csv(file, index=False)

    return path
