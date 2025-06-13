import numpy as np
import pandas as pd
import xarray as xr
from astropy import units as u


def calculate_mass_fractions_over_time(isotope_mass_fraction, time_coord):
    """
    Calculate the evolution of isotope mass fractions over time due to radioactive decay.

    Uses the existing calculate_decayed_mass_fractions method on the isotope_mass_fraction
    object to compute decay for each time step, then reassembles as a DataArray.

    Parameters
    ----------
    isotope_mass_fraction : DataFrame-like object
        Object with isotope mass fractions that has a calculate_decayed_mass_fractions method.
        Index should be MultiIndex with (atomic_number, mass_number).
    time_coord : xr.DataArray
        Time coordinate dataset containing time values.

    Returns
    -------
    xr.DataArray
        DataArray with dimensions (isotope, cell_id, time_step_id) containing evolved mass fractions.
        Includes atomic_number and mass_number coordinates along the isotope dimension.
    """
    # Get time step information from the time coordinate
    time_step_ids = time_coord.time_step_id.values
    delta_times = time_coord.delta_time.values * u.Unit(time_coord.attrs["unit"])

    # Get cell IDs from the original object
    cell_ids = isotope_mass_fraction.columns

    # Initialize with the undecayed mass fractions
    current_mass_fractions = isotope_mass_fraction.copy()
    all_isotopes = set(current_mass_fractions.index)
    mass_fractions_list = [current_mass_fractions]  # Start with undecayed

    # Iteratively decay for each delta_time (ignore the last entry)
    for delta_time in delta_times[:-1]:
        current_mass_fractions = (
            current_mass_fractions.calculate_decayed_mass_fractions(delta_time)
        )
        mass_fractions_list.append(current_mass_fractions)
        all_isotopes.update(current_mass_fractions.index)

    # Create unified isotope index from all isotopes found across all time steps
    isotope_index = pd.MultiIndex.from_tuples(
        sorted(all_isotopes), names=["atomic_number", "mass_number"]
    )

    # Second pass: reindex all time steps to the unified isotope index and collect values
    reindexed_mass_fractions_list = []
    for mass_fractions in mass_fractions_list:
        # Reindex to include all isotopes, filling missing ones with 0
        reindexed_fractions = mass_fractions.reindex(isotope_index, fill_value=0.0)
        reindexed_mass_fractions_list.append(reindexed_fractions.values)

    # Stack the results along the time dimension
    mass_fractions_array = np.stack(reindexed_mass_fractions_list, axis=-1)

    # Create output DataArray with proper coordinates
    mass_fractions_over_time = xr.DataArray(
        mass_fractions_array,
        dims=["isotope", "cell_id", "time_step_id"],
        coords={
            "isotope": isotope_index,
            "cell_id": cell_ids,
            "time_step_id": time_step_ids,
            "time_step_start": ("time_step_id", time_coord.time_step_start.values),
            "time_step_stop": ("time_step_id", time_coord.time_step_stop.values),
        },
    )

    return mass_fractions_over_time


def calculate_decays_over_time(isotope_mass_fraction, time_coord, cell_masses=None):
    """
    Calculate the number of decays for each isotope over time due to radioactive decay for each time_step.

    Uses the existing calculate_number_of_decays method on the isotope_mass_fraction
    object to compute total decays for each time step, then reassembles as a DataArray.

    Parameters
    ----------
    isotope_mass_fraction : DataFrame-like object
        Object with isotope mass fractions that has a calculate_number_of_decays method.
        Index should be MultiIndex with (atomic_number, mass_number).
    time_coord : xr.DataArray
        Time coordinate dataset containing time values.
    cell_masses : astropy.units.Quantity, optional
        1D array of shell masses. If provided, isotopic abundances will be
        multiplied by the corresponding shell mass.

    Returns
    -------
    xr.DataArray
        DataArray with dimensions (isotope, cell_id, time_step_id) containing total decays.
        Includes atomic_number and mass_number coordinates along the isotope dimension.
    """
    # Get time step information from the time coordinate
    time_step_ids = time_coord.time_step_id.values
    delta_times = time_coord.delta_time.values * u.Unit(time_coord.attrs["unit"])

    # Get cell IDs from the original object
    cell_ids = isotope_mass_fraction.columns

    # Initialize with the undecayed mass fractions
    current_mass_fractions = isotope_mass_fraction.copy()
    all_isotopes = set(current_mass_fractions.index)
    total_decays_list = []

    # Calculate total decays for each delta_time (ignore the last entry)
    for delta_time in delta_times:
        # Calculate total decays for this time step
        total_decays = current_mass_fractions.calculate_number_of_decays(
            delta_time, cell_masses
        )
        total_decays_list.append(total_decays)
        all_isotopes.update(total_decays.index)

        # Update mass fractions for next iteration
        current_mass_fractions = (
            current_mass_fractions.calculate_decayed_mass_fractions(delta_time)
        )

    # Create unified isotope index from all isotopes found across all time steps
    isotope_index = pd.MultiIndex.from_tuples(
        sorted(all_isotopes), names=["atomic_number", "mass_number"]
    )

    # Second pass: reindex all time steps to the unified isotope index and collect values
    reindexed_decays_list = []
    for total_decays in total_decays_list:
        # Reindex to include all isotopes, filling missing ones with 0
        reindexed_decays = total_decays.reindex(isotope_index, fill_value=0.0)
        reindexed_decays_list.append(reindexed_decays.values)

    # Stack the results along the time dimension
    total_decays_array = np.stack(reindexed_decays_list, axis=-1)

    # Create output DataArray with proper coordinates
    total_decays_over_time = xr.DataArray(
        total_decays_array,
        dims=["isotope", "cell_id", "time_step_id"],
        coords={
            "isotope": isotope_index,
            "cell_id": cell_ids,
            "time_step_id": time_step_ids,
            "time_step_start": ("time_step_id", time_coord.time_step_start.values),
            "time_step_stop": ("time_step_id", time_coord.time_step_stop.values),
        },
    )

    return total_decays_over_time
