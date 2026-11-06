import pandas as pd
from astropy import units as u

from tardis.io.hdf_writer_mixin import HDFWriterMixin
from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.transport.montecarlo.packets.packet_collections import (
    PacketCollection,
)


class MonteCarloTransportState(HDFWriterMixin):
    hdf_properties = [
        "output_nu",
        "output_energy",
        "nu_bar_estimator",
        "j_estimator",
        "j_blue_estimator",
        "packet_luminosity",
        "time_of_simulation",
        "emitted_packet_mask",
        "last_interaction_type",
        "last_interaction_in_nu",
        "last_interaction_in_r",
        "last_line_interaction_out_id",
        "last_line_interaction_in_id",
        "last_line_interaction_shell_id",
    ]

    hdf_name = "transport_state"

    last_interaction_type = None
    last_interaction_in_nu = None
    last_interaction_in_r = None
    last_line_interaction_out_id = None
    last_line_interaction_in_id = None
    last_line_interaction_shell_id = None

    def __init__(
        self,
        packet_collection: PacketCollection,
        geometry_state_numba: NumbaRadial1DGeometry,
        opacity_state_numba: OpacityStateNumba,
        time_explosion: u.Quantity,
        n_levels_bf_species_by_n_cells_tuple: tuple[int, int],
        tracker_full_df: pd.DataFrame | None = None,
        tracker_last_interaction_df: pd.DataFrame | None = None,
    ) -> None:
        """
        Initialize the state used by Monte Carlo radiative transfer.

        Parameters
        ----------
        packet_collection : tardis.transport.montecarlo.packets.packet_collections.PacketCollection
            Collection of packets propagated during the simulation.
        geometry_state_numba : tardis.model.geometry.radial1d.NumbaRadial1DGeometry
            Numba-compatible geometry state describing the simulation shells.
        opacity_state_numba : tardis.opacities.opacity_state_numba.OpacityStateNumba
            Numba-compatible opacity state used for packet interactions.
        time_explosion : astropy.units.Quantity
            Time since the explosion in seconds.
        n_levels_bf_species_by_n_cells_tuple : tuple[int, int]
            Shape of the bound-free level and cell dimensions.
        tracker_full_df : pandas.DataFrame, optional
            DataFrame containing full packet interaction tracking data.
        tracker_last_interaction_df : pandas.DataFrame, optional
            DataFrame containing the last interaction for each packet.
        """
        self.packet_collection = packet_collection
        self.n_levels_bf_species_by_n_cells_tuple = (
            n_levels_bf_species_by_n_cells_tuple
        )
        self.estimators_bulk = None
        self.estimators_line = None
        self.estimators_continuum = None
        self.enable_full_relativity = False
        self.enable_continuum_processes = False
        self.time_explosion = time_explosion
        self.geometry_state_numba = geometry_state_numba
        self.opacity_state_numba = opacity_state_numba
        self.tracker_full_df = tracker_full_df
        self.tracker_last_interaction_df = tracker_last_interaction_df

    @property
    def output_nu(self):
        return self.packet_collection.output_nus * u.Hz

    @property
    def output_energy(self):
        return self.packet_collection.output_energies * u.erg

    @property
    def nu_bar_estimator(self):
        return self.estimators_bulk.mean_frequency

    @property
    def j_estimator(self):
        return self.estimators_bulk.mean_intensity_total

    @property
    def j_blue_estimator(self):
        return self.estimators_line.mean_intensity_blueward

    @property
    def time_of_simulation(self):
        return self.packet_collection.time_of_simulation * u.s

    @property
    def packet_luminosity(self):
        return (
            self.packet_collection.output_energies
            * u.erg
            / (self.packet_collection.time_of_simulation * u.s)
        )

    @property
    def emitted_packet_mask(self):
        return self.packet_collection.output_energies >= 0

    @property
    def emitted_packet_nu(self):
        return (
            self.packet_collection.output_nus[self.emitted_packet_mask] * u.Hz
        )

    @property
    def reabsorbed_packet_nu(self):
        return (
            self.packet_collection.output_nus[~self.emitted_packet_mask] * u.Hz
        )

    @property
    def emitted_packet_luminosity(self):
        return self.packet_luminosity[self.emitted_packet_mask]

    @property
    def reabsorbed_packet_luminosity(self):
        return -self.packet_luminosity[~self.emitted_packet_mask]
