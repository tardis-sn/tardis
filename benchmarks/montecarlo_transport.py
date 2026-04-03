"""
Basic TARDIS Benchmark.
"""

from numba.typed import List

from benchmarks.benchmark_base import BenchmarkBase
from tardis.transport.montecarlo.modes.classic.montecarlo_transport import (
    montecarlo_transport,
)
from tardis.transport.montecarlo.packets.packet_collections import (
    PacketCollection,
)
from tardis.transport.montecarlo.packets.trackers.tracker_last_interaction import (
    TrackerLastInteraction as _TrackerLastInteraction,
)


class BenchmarkMonteCarloTransport(BenchmarkBase):
    """
    Benchmarks for the classic Monte Carlo transport workflow.
    """

    repeat = 2

    def setup(self):
        self.nb_simulation = self.nb_simulation_verysimple

        transport_state = self.nb_simulation.transport.transport_state
        simulation_state = self.nb_simulation.simulation_state
        plasma_state = self.nb_simulation.plasma
        packet_collection = transport_state.packet_collection

        self.initial_radii = packet_collection.initial_radii.copy()
        self.initial_nus = packet_collection.initial_nus.copy()
        self.initial_mus = packet_collection.initial_mus.copy()
        self.initial_energies = packet_collection.initial_energies.copy()
        self.packet_seeds = packet_collection.packet_seeds.copy()
        self.radiation_field_luminosity = packet_collection.radiation_field_luminosity
        self.no_of_packets = len(self.initial_nus)

        geometry_state = getattr(transport_state, "geometry_state", None)
        if geometry_state is None:
            geometry_state = simulation_state.geometry
            if hasattr(geometry_state, "to_numba"):
                geometry_state = geometry_state.to_numba()
        self.geometry_state = geometry_state

        time_explosion = getattr(transport_state, "time_explosion", None)
        if time_explosion is None:
            time_explosion = simulation_state.time_explosion
        self.time_explosion = time_explosion.to("s").value

        opacity_state = getattr(transport_state, "opacity_state", None)
        if opacity_state is None:
            opacity_state = getattr(plasma_state, "opacity_state", None)
            if opacity_state is None:
                raise AttributeError(
                    "Could not find opacity state on transport_state or plasma_state"
                )
            if hasattr(opacity_state, "to_numba"):
                opacity_state = opacity_state.to_numba()
        self.opacity_state = opacity_state

        self.montecarlo_configuration = (
            self.nb_simulation.transport.montecarlo_configuration
        )
        self.spectrum_frequency_grid = (
            self.nb_simulation.transport.spectrum_frequency_grid.value.copy()
        )
        self.number_of_vpackets = self.montecarlo_configuration.NUMBER_OF_VPACKETS
        self.show_progress_bars = False

    def _make_packet_collection(self):
        return PacketCollection(
            self.initial_radii.copy(),
            self.initial_nus.copy(),
            self.initial_mus.copy(),
            self.initial_energies.copy(),
            self.packet_seeds.copy(),
            self.radiation_field_luminosity,
        )

    def _make_trackers(self, no_of_packets):
        trackers = List()
        for _ in range(no_of_packets):
            trackers.append(_TrackerLastInteraction())
        return trackers

    def time_montecarlo_transport(self):
        packet_collection = self._make_packet_collection()
        trackers = self._make_trackers(self.no_of_packets)

        montecarlo_transport(
            packet_collection,
            self.geometry_state,
            self.time_explosion,
            self.opacity_state,
            self.montecarlo_configuration,
            self.spectrum_frequency_grid,
            trackers,
            self.number_of_vpackets,
            self.show_progress_bars,
        )
