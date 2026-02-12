import warnings

from astropy import units as u

from tardis.io.hdf_writer_mixin import HDFWriterMixin


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

    vpacket_hdf_properties = [
        "virt_packet_nus",
        "virt_packet_energies",
        "virt_packet_initial_rs",
        "virt_packet_initial_mus",
        "virt_packet_last_interaction_in_nu",
        "virt_packet_last_interaction_in_r",
        "virt_packet_last_interaction_type",
        "virt_packet_last_line_interaction_in_id",
        "virt_packet_last_line_interaction_out_id",
        "virt_packet_last_line_interaction_shell_id",
    ]

    hdf_name = "transport_state"

    last_interaction_type = None
    last_interaction_in_nu = None
    last_interaction_in_r = None
    last_line_interaction_out_id = None
    last_line_interaction_in_id = None
    last_line_interaction_shell_id = None

    virt_logging = False

    def __init__(
        self,
        packet_collection,
        geometry_state,
        opacity_state,
        time_explosion,
        n_levels_bf_species_by_n_cells_tuple,
        tracker_full_df=None,
        tracker_last_interaction_df=None,
        vpacket_tracker=None,
    ):
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
        self.geometry_state = geometry_state
        self.opacity_state = opacity_state
        self.tracker_full_df = tracker_full_df
        self.tracker_last_interaction_df = tracker_last_interaction_df
        self.vpacket_tracker = vpacket_tracker

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

    @property
    def virt_packet_nus(self):
        try:
            return u.Quantity(self.vpacket_tracker.nus, u.Hz)
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virt_packet_nus:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virt_packet_energies(self):
        try:
            return u.Quantity(self.vpacket_tracker.energies, u.erg)
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virt_packet_energies:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virtual_packet_luminosity(self):
        try:
            return (
                self.virt_packet_energies
                / self.packet_collection.time_of_simulation
            )
        except TypeError:
            warnings.warn(
                "MontecarloTransport.virtual_packet_luminosity:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virt_packet_initial_rs(self):
        try:
            return u.Quantity(self.vpacket_tracker.initial_rs, u.erg)
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virt_packet_initial_rs:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virt_packet_initial_mus(self):
        try:
            return u.Quantity(self.vpacket_tracker.initial_mus, u.erg)
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virt_packet_initial_mus:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virt_packet_last_interaction_in_nu(self):
        try:
            return u.Quantity(self.vpacket_tracker.last_interaction_in_nu, u.Hz)
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virt_packet_last_interaction_in_nu:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virt_packet_last_interaction_in_r(self):
        try:
            return u.Quantity(self.vpacket_tracker.last_interaction_in_r, u.cm)
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virt_packet_last_interaction_in_r:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virt_packet_last_interaction_type(self):
        try:
            return self.vpacket_tracker.last_interaction_type
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virt_packet_last_interaction_type:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virt_packet_last_line_interaction_in_id(self):
        try:
            return self.vpacket_tracker.last_interaction_in_id
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virt_packet_last_line_interaction_in_id:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virt_packet_last_line_interaction_out_id(self):
        try:
            return self.vpacket_tracker.last_interaction_out_id
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virt_packet_last_line_interaction_out_id:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virt_packet_last_line_interaction_shell_id(self):
        try:
            return self.vpacket_tracker.last_interaction_shell_id
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virt_packet_last_line_interaction_shell_id:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None
