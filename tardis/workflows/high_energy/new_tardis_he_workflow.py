import logging
from dataclasses import dataclass

import astropy.units as u
import numpy as np
import xarray as xr

from tardis.energy_input.decay_over_time import calculate_decays_over_time
from tardis.energy_input.decay_radiation import (
    get_decay_radiation_data,
    convert_decay_radiation_to_ds,
)
from tardis.energy_input.energy_source import get_radioactive_isotopes
from tardis.energy_input.util import create_time_coord_from_array
from tardis.model import SimulationState

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class TARDISHEWorkflow:
    def __init__(self, atom_data, configuration, config_type="yaml"):

        if config_type == "csvy":
            self.simulation_state = SimulationState.from_csvy(configuration, atom_data)
        else:
            self.simulation_state = SimulationState.from_config(
                configuration, atom_data
            )

        self.shell_masses = self.simulation_state.volume * self.simulation_state.density

        self.isotopic_mass_fraction = (
            self.simulation_state.composition.isotopic_mass_fraction
        )

        # Load and convert decay radiation data (already filtered to radioactive isotopes)
        em_radiation_data, bp_radiation_data = get_decay_radiation_data(
            atom_data.decay_radiation_data, self.isotopic_mass_fraction.index
        )

        # Convert to xarray datasets and save to self
        self.em_radiation_ds = convert_decay_radiation_to_ds(em_radiation_data)
        self.bp_radiation_ds = convert_decay_radiation_to_ds(bp_radiation_data)

        self.radioactive_isotopes = get_radioactive_isotopes(
            self.isotopic_mass_fraction.index
        )

    def run(
        self,
        time_start,
        time_end,
        time_steps,
        time_space="linear",
    ):
        """
        Run the workflow to calculate total decays per isotope over time.

        Parameters
        ----------
        time_start : float
            The start time in days.
        time_end : float
            The end time in days.
        time_steps : int
            The number of time steps to take between time_start and time_end.
        time_space : str, optional
            The time spacing between time steps ('linear' or 'log'), by default 'linear'

        Returns
        -------
        TARDISHEWorkflowResult
            Result containing total decays over time
        """
        # Create time array
        if time_space == "log":
            time_values = np.geomspace(time_start, time_end, time_steps + 1) * u.day
        else:
            time_values = np.linspace(time_start, time_end, time_steps + 1) * u.day

        # Create time coordinate
        time_coord = create_time_coord_from_array(time_values, unit=u.day)

        decays_over_time = calculate_decays_over_time(
            self.isotopic_mass_fraction, time_coord, cell_masses=self.shell_masses
        )

        # Calculate total decay energy by aligning and multiplying datasets
        decay_energy_over_time = self.calculate_decay_energy(decays_over_time)

        return TARDISHEWorkflowResult(
            decays_over_time=decays_over_time,
            decay_energy_over_time=decay_energy_over_time,
        )

    def calculate_decay_energy(self, decays_over_time):
        """
        Calculate total decay energy by multiplying em_radiation data with number of decays.

        Parameters
        ----------
        decays_over_time : xr.DataArray
            The total number of decays over time for each isotope and cell

        Returns
        -------
        xr.DataArray
            Total decay energy in keV for each decay radiation channel, cell, and time step
        """
        # Get the multi-index for isotopes from em_radiation_ds
        multi_index = self.em_radiation_ds.indexes["isotope"]

        # Select decays for isotopes that have radiation data
        broadcast_decays_over_time = decays_over_time.sel(isotope=multi_index)

        # Remove isotope-related coordinates and set decay_radiation_channel_id as coordinate
        broadcast_decays_over_time = broadcast_decays_over_time.drop_vars(
            ["isotope", "atomic_number", "mass_number"]
        ).reset_coords(drop=True)

        # Add decay_radiation_channel_id as a coordinate
        broadcast_decays_over_time = broadcast_decays_over_time.assign_coords(
            isotope=self.em_radiation_ds.decay_radiation_channel_id.values
        )

        # Set decay_radiation_channel_id as an index for efficient selection
        broadcast_decays_over_time = broadcast_decays_over_time.rename(
            {"isotope": "decay_radiation_channel_id"}
        )

        decay_energy_over_time = (
            broadcast_decays_over_time * self.em_radiation_ds.energy_per_decay_kev
        )
        decay_energy_over_time = decay_energy_over_time.set_xindex(
            ["atomic_number", "mass_number"]
        )

        decay_energy_over_time.name = "total_decay_energy_kev"
        decay_energy_over_time.attrs["unit"] = "keV"
        return decay_energy_over_time


@dataclass
class TARDISHEWorkflowResult:
    decays_over_time: xr.DataArray
    decay_energy_over_time: xr.DataArray
