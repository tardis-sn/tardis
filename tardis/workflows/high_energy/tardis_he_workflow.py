import logging
from dataclasses import dataclass

import pandas as pd
from astropy import units as u

from tardis.energy_input.decay_radiation import get_decay_radiation_data
from tardis.energy_input.gamma_ray_channel import (
    calculate_total_decays,
    create_inventories_dict,
    create_isotope_decay_df,
    create_isotope_dicts,
    time_evolve_cumulative_decay,
)
from tardis.energy_input.main_gamma_ray_loop import (
    get_effective_time_array,
    run_gamma_ray_loop,
)
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

        self.atom_data = atom_data
        self.gamma_ray_lines = atom_data.decay_radiation_data

        self.isotopic_mass_fraction = (
            self.simulation_state.composition.isotopic_mass_fraction
        )

    def calculate_total_decays(self, time_start, time_end):
        """
        Calculate the total number of decays for each isotope in the simulation
        state between t_start and t_end for each shell.

        Parameters
        ----------
        time_start : float
            The start time in days.
        time_end : float
            The end time in days.

        Returns
        -------
        total_decays : dataframe
            A dataframe containing the total number of decays for each isotope
            in the simulation state between t_start and t_end for each shell.
        """
        isotopes = create_isotope_dicts(self.isotopic_mass_fraction, self.cell_masses)
        inventories = create_inventories_dict(isotopes)

        total_decays = calculate_total_decays(inventories, time_end - time_start)

        return total_decays

    def decay_isotopes_expanded(self, total_decays):
        """
        Expand the total decays dataframe to include the decay energy for each gamma-X ray
        transition in each shell between t_start and t_end.

        Parameters
        ----------
        total_decays : dataframe
            A dataframe containing the total number of decays for each isotope
            in the simulation state between t_start and t_end for each shell.

        Returns
        -------
        decay_isotopes : dataframe
            A dataframe containing the total number of decays for each isotope
            in the simulation state between t_start and t_end for each shell
            expanded to include the decay energy for each gamma-X ray transition.
        """
        return create_isotope_decay_df(total_decays, self.gamma_ray_lines)

    def time_evolve_cumulative_decay_expanded(self, times):
        """
        Time evolve the cumulative decay for each isotope in the simulation state
        between t_start and t_end for each shell.

        Parameters
        ----------
        times : array
            An array of times in days.

        Returns
        -------
        cumulative_decay : dataframe
            A dataframe containing the cumulative decay for each isotope, each transition
            each shell for all time steps between t_start and t_end.
        """
        return time_evolve_cumulative_decay(
            self.isotopic_mass_fraction, self.cell_masses, self.gamma_ray_lines, times
        )

    def get_times(self, time_start, time_end, time_space, time_steps):
        """
        Create an array of times between t_start and t_end.

        Parameters
        ----------
        time_start : float
            The start time in days.
        time_end : float
            The end time in days.
        time_space : string
            The time spacing between time steps ('linear' or 'log').
        time_steps : int
            The number of time steps to take between t_start and t_end.

        Returns
        -------
        times : array
            An array of times in days.
        """
        return get_effective_time_array(time_start, time_end, time_space, time_steps)

    def run(
        self,
        time_start,
        time_end,
        number_of_packets,
        time_steps,
        time_space,
        seed,
        fp,
        spectrum_bins,
        grey_opacity=-1,
        legacy=False,
        legacy_atom_data=None,
    ):
        """
        Run the gamma-ray transport simulation.

        Parameters
        ----------
        time_start : float
            The start time in days.
        time_end : float
            The end time in days.
        times : array
            An array of times in days.
        number_of_packets : int
            The number of packets to simulate.
        time_steps : int
            The number of time steps to take between t_start and t_end.
        time_space : string
            The time spacing between time steps.
        seed : int
            The random seed for the simulation.
        fp : float
            Positronium fraction
        spectrum_bins : int
            The number of bins in the spectrum.
        grey_opacity : float
            The grey opacity of the simulation.
        """
        time_start = u.Quantity(time_start, u.day)
        time_end = u.Quantity(time_end, u.day)

        times, effective_times = self.get_times(
            time_start.to(u.day).value, time_end.to(u.day).value, time_space, time_steps
        )

        cell_masses = (self.simulation_state.volume * self.simulation_state.density).to(
            u.g
        )
        self.cell_masses = cell_masses  # remove once the total decays works
        total_decays = self.calculate_total_decays(time_start, time_end)

        total_decays2 = self.simulation_state.composition.isotopic_mass_fraction.calculate_number_of_decays(
            time_end - time_start, cell_masses
        )
        decay_isotopes = self.decay_isotopes_expanded(total_decays)

        em_radiation_data, bp_radiation_data = get_decay_radiation_data(
            self.atom_data.decay_radiation_data,
            self.simulation_state.composition.isotopic_mass_fraction.index,
        )
        decay_over_time = self.time_evolve_cumulative_decay_expanded(times)

        (
            escape_energy,
            escape_energy_cosi,
            packets_escaped,
            gamma_ray_deposited_energy,
            total_deposited_energy,
            positron_energy,
        ) = run_gamma_ray_loop(
            self.simulation_state,
            decay_isotopes,
            decay_over_time,
            number_of_packets,
            times,
            effective_times,
            seed,
            fp,
            spectrum_bins,
            grey_opacity,
            legacy=legacy,
            legacy_atom_data=legacy_atom_data,
        )

        return TARDISHEWorkflowResult(
            escape_energy=escape_energy,
            escape_energy_cosi=escape_energy_cosi,
            packets_escaped=packets_escaped,
            gamma_ray_deposited_energy=gamma_ray_deposited_energy,
            total_deposited_energy=total_deposited_energy,
            positron_energy=positron_energy,
        )


@dataclass
class TARDISHEWorkflowResult:

    escape_energy: pd.DataFrame
    escape_energy_cosi: pd.DataFrame
    packets_escaped: pd.DataFrame
    gamma_ray_deposited_energy: pd.DataFrame
    total_deposited_energy: pd.DataFrame
    positron_energy: pd.DataFrame
