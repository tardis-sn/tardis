from __future__ import annotations

from dataclasses import dataclass

import astropy.units as u
from pathlib import Path
from tardis.visualization import plot_util as pu


@dataclass
class SDECState:
    """
    Container for all data required to perform SDEC analysis and plotting.
    """

    packet_data: dict[str, object]
    time_of_simulation: u.Quantity
    plot_wavelength: u.Quantity
    plot_frequency: u.Quantity
    plot_frequency_bins: u.Quantity
    modeled_spectrum_luminosity: u.Quantity
    spectrum_data: dict[str, u.Quantity]

    @classmethod
    def from_simulation(
        cls, simulation, packets_mode: str = "real", packet_wvl_range=None
    ) -> SDECState:
        if packet_wvl_range is None:
            packet_wvl_range = [3000, 9000] * u.AA
        packet_data = pu.extract_and_process_packet_data(
            simulation.transport.transport_state,
            simulation.plasma,
            packets_mode,
        )
        time_of_simulation = (
            simulation.transport.transport_state.packet_collection.time_of_simulation
            * u.s
        )
        plot_frequency_bins = getattr(
            simulation.spectrum_solver, f"spectrum_{packets_mode}_packets"
        )._frequency
        plot_wavelength = getattr(
            simulation.spectrum_solver, f"spectrum_{packets_mode}_packets"
        ).wavelength
        modeled_spectrum_luminosity = getattr(
            simulation.spectrum_solver, f"spectrum_{packets_mode}_packets"
        ).luminosity_density_lambda
        plot_frequency = plot_frequency_bins[:-1]

        spectrum_data = pu.get_spectrum_data_from_spectrum_solver(
            simulation.spectrum_solver, packets_mode=packets_mode
        )
        return cls(
            packet_data=packet_data,
            time_of_simulation=time_of_simulation,
            plot_wavelength=plot_wavelength,
            plot_frequency=plot_frequency,
            plot_frequency_bins=plot_frequency_bins,
            modeled_spectrum_luminosity=modeled_spectrum_luminosity,
            spectrum_data=spectrum_data,
        )

    @classmethod
    def from_workflow(
        cls, workflow, packets_mode: str = "real", packet_wvl_range=None
    ) -> SDECState:
        if packet_wvl_range is None:
            packet_wvl_range = [3000, 9000] * u.AA
        packet_data = pu.extract_and_process_packet_data(
            workflow.transport_state,
            workflow.plasma_solver,
            packets_mode,
        )
        time_of_simulation = (
            workflow.transport_state.packet_collection.time_of_simulation * u.s
        )
        plot_frequency_bins = getattr(
            workflow.spectrum_solver, f"spectrum_{packets_mode}_packets"
        )._frequency
        plot_wavelength = getattr(
            workflow.spectrum_solver, f"spectrum_{packets_mode}_packets"
        ).wavelength
        modeled_spectrum_luminosity = getattr(
            workflow.spectrum_solver, f"spectrum_{packets_mode}_packets"
        ).luminosity_density_lambda
        plot_frequency = plot_frequency_bins[:-1]
        spectrum_data = pu.get_spectrum_data_from_spectrum_solver(
            workflow.spectrum_solver, packets_mode=packets_mode
        )
        return cls(
            packet_data=packet_data,
            time_of_simulation=time_of_simulation,
            plot_wavelength=plot_wavelength,
            plot_frequency=plot_frequency,
            plot_frequency_bins=plot_frequency_bins,
            modeled_spectrum_luminosity=modeled_spectrum_luminosity,
            spectrum_data=spectrum_data,
        )

    @classmethod
    def from_workflow_hdf(
        cls,
        hdf_store,
        lines_hdf_path,
        transport_state_hdf_path,
        spectrum_solver_hdf_path,
        packets_mode: str = "real",
    ) -> SDECState:
        if isinstance(hdf_store, str):
            import pandas as pd

            hdf_store = pd.HDFStore(hdf_store, "r")
            close_hdf = True
        else:
            close_hdf = False
        try:
            # Use plot_util to extract packet data
            packet_data = pu.extract_and_process_packet_data_hdf(
                hdf_store,
                packets_mode,
                lines_hdf_path=lines_hdf_path,
                transport_state_hdf_path=transport_state_hdf_path,
            )
            # Extract time_of_simulation
            time_of_simulation = (
                hdf_store[str(Path("/", transport_state_hdf_path, "scalars"))][
                    "time_of_simulation"
                ][0]
                * u.s
            )
            spectrum_data = pu.get_spectrum_data_from_hdf(
                hdf_store,
                packets_mode=packets_mode,
                spectrum_solver_hdf_path=spectrum_solver_hdf_path,
            )
        finally:
            if close_hdf:
                hdf_store.close()
        return cls(
            packet_data=packet_data,
            time_of_simulation=time_of_simulation,
            plot_wavelength=plot_wavelength,
            plot_frequency=plot_frequency,
            plot_frequency_bins=plot_frequency_bins,
            modeled_spectrum_luminosity=modeled_spectrum_luminosity,
            spectrum_data=spectrum_data,
        )
