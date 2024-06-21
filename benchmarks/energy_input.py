import astropy.units as u
import numpy as np

from benchmarks.benchmark_base import BenchmarkBase
from tardis.energy_input import GXPacket, main_gamma_ray_loop
from tardis.energy_input.energy_source import (
    get_nuclear_lines_database,
)
from tardis.energy_input.gamma_packet_loop import gamma_packet_loop
from tardis.energy_input.gamma_ray_packet_source import RadioactivePacketSource
from tardis.energy_input.gamma_ray_transport import (
    calculate_average_energies,
    calculate_average_power_per_mass,
    calculate_ejecta_velocity_volume,
    calculate_energy_per_mass,
    calculate_total_decays_old,
    create_inventories_dict_old,
    create_isotope_dicts_old,
    decay_chain_energies,
    distribute_packets,
    get_taus,
    iron_group_fraction_per_shell,
)
from tardis.io.atom_data import AtomData
from tardis.io.configuration import config_reader
from tardis.model import SimulationState
from tardis.plasma.base import BasePlasma
from tardis.plasma.properties import (
    Abundance,
    AtomicData,
    AtomicMass,
    Density,
    IsotopeAbundance,
    IsotopeMass,
    IsotopeNumberDensity,
    NumberDensity,
    SelectedAtoms,
)


class BenchmarkMainGammaRayLoop(BenchmarkBase):
    """
    Class to benchmark the gamma packet loop function.
    """

    def setup(self):
        self.atom_data_file = "tardis/io/configuration/tests/data/tardis_configv1_density_exponential_nebular_Ni_only.yml"
        atom_data = AtomData.from_hdf(self.atom_data_file)

        config = config_reader.Configuration.from_yaml(self.atom_data_file)

        self.model = SimulationState.from_config(config, atom_data)

        input = [Density, Abundance, IsotopeAbundance, AtomicData, AtomicMass, IsotopeNumberDensity, NumberDensity, SelectedAtoms, IsotopeMass]
        self.plasma = BasePlasma(
            plasma_properties=input,
            density = self.model.density,
            abundance=self.model.abundance,
            isotope_abundance=self.model.composition.raw_isotope_abundance,
            atomic_data = atom_data
        )

        self.num_packets =  100000
        self.time_start = 0.0011574074
        self.time_end = 20.0
        self.time_space = "log"
        self.time_steps = 50
        self.seed = 1
        self.positronium_fraction = 0.0
        self.spectrum_bins = 1000
        self.grey_opacity = -1
        self.photoabsorption_opacity = "tardis"
        self.pair_creation_opacity = "tardis"


    def time_main_gamma_ray_loop(self):
        main_gamma_ray_loop.run_gamma_ray_loop(
            self.model,
            self.plasma,
            num_decays=self.num_packets,
            time_start=0.0011574074,
            time_end=20.0,
            time_space="log",
            time_steps=50,
            seed=1,
            positronium_fraction=0.0,
            spectrum_bins=1000,
            grey_opacity=-1,
            path_to_decay_data=self.atom_data_file
        )



class BenchmarkGammaPacketLoop(BenchmarkMainGammaRayLoop):
    """
    Class to benchmark the gamma packet loop function.
    """

    def setup(self):
        super().setup()
        np.random.seed(seed=1)
        self.inner_velocities = self.model.v_inner.to("cm/s").value
        self.outer_velocities = self.model.v_outer.to("cm/s").value
        ejecta_volume = self.model.volume.to("cm^3").value
        number_of_shells = self.model.no_of_shells
        shell_masses = self.model.volume * self.model.density
        raw_isotope_abundance = self.model.composition.raw_isotope_abundance.sort_values(
            by=["atomic_number", "mass_number"], ascending=False
        )
        self.time_start *= u.d.to(u.s)
        self.time_end *= u.d.to(u.s)

        if self.time_space == "log":
            self.times = np.geomspace(self.time_start, self.time_end, self.time_steps + 1)
            self.effective_time_array = np.sqrt(self.times[:-1] * self.times[1:])
        else:
            self.times = np.linspace(self.time_start, self.time_end, self.time_steps + 1)
            self.effective_time_array = 0.5 * (self.times[:-1] + self.times[1:])

        self.dt_array = np.diff(self.times)

        ejecta_velocity_volume = calculate_ejecta_velocity_volume(self.model)

        self.inv_volume_time = (
            1.0 / ejecta_velocity_volume[:, np.newaxis]
        ) / self.effective_time_array**3.0

        self.energy_df_rows = np.zeros((number_of_shells, self.time_steps))
        for atom_number in self.plasma.isotope_number_density.index.get_level_values(0):
            values = self.plasma.isotope_number_density.loc[atom_number].values
            if values.shape[0] > 1:
                self.plasma.isotope_number_density.loc[atom_number].update = np.sum(
                    values, axis=0
                )
            else:
                self.plasma.isotope_number_density.loc[atom_number].update = values

        electron_number_density = self.plasma.number_density.mul(
            self.plasma.number_density.index, axis=0
        ).sum()
        taus, parents = get_taus(raw_isotope_abundance)
        electron_number = np.array(electron_number_density * ejecta_volume)
        self.electron_number_density_time = (
            electron_number[:, np.newaxis] * self.inv_volume_time
        )

        self.mass_density_time = shell_masses[:, np.newaxis] * self.inv_volume_time
        gamma_ray_lines = get_nuclear_lines_database(self.atom_data_file)
        isotope_dict = create_isotope_dicts_old(raw_isotope_abundance, shell_masses)
        inventories_dict = create_inventories_dict_old(isotope_dict)
        total_decays = calculate_total_decays_old(
            inventories_dict, self.time_end - self.time_start
        )

        (
            average_energies,
            average_positron_energies,
            gamma_ray_line_dict,
        ) = calculate_average_energies(raw_isotope_abundance, gamma_ray_lines)

        decayed_energy = decay_chain_energies(
            average_energies,
            total_decays,
        )
        energy_per_mass, energy_df = calculate_energy_per_mass(
            decayed_energy, raw_isotope_abundance, shell_masses
        )
        average_power_per_mass = calculate_average_power_per_mass(
            energy_per_mass, self.time_end - self.time_start
        )

        total_energy = energy_df.sum().sum()
        packets_per_isotope_df = (
            distribute_packets(decayed_energy, total_energy, self.num_packets)
            .round()
            .fillna(0)
            .astype(int)
        )

        total_energy = total_energy * u.eV.to("erg")

        self.iron_group_fraction_per_shell = iron_group_fraction_per_shell(self.model).to_numpy()
        number_of_packets = packets_per_isotope_df.sum().sum()
        individual_packet_energy = total_energy / number_of_packets

        packet_source = RadioactivePacketSource(
            individual_packet_energy,
            gamma_ray_line_dict,
            self.positronium_fraction,
            self.inner_velocities,
            self.outer_velocities,
            self.inv_volume_time,
            self.times,
            self.energy_df_rows,
            self.effective_time_array,
            taus,
            parents,
            average_positron_energies,
            average_power_per_mass,
        )

        packet_collection = packet_source.create_packets(packets_per_isotope_df)

        self.energy_df_rows = packet_source.energy_df_rows
        self.energy_plot_df_rows = np.zeros((number_of_packets, 8))

        self.packets = []
        for i in range(number_of_packets):
            packet = GXPacket(
                packet_collection.location[:, i],
                packet_collection.direction[:, i],
                packet_collection.energy_rf[i],
                packet_collection.energy_cmf[i],
                packet_collection.nu_rf[i],
                packet_collection.nu_cmf[i],
                packet_collection.status[i],
                packet_collection.shell[i],
                packet_collection.time_current[i],
            )
            self.packets.append(packet)
            self.energy_plot_df_rows[i] = np.array(
                [
                    i,
                    packet.energy_rf,
                    packet.get_location_r(),
                    packet.time_current,
                    int(packet.status),
                    0,
                    0,
                    0,
                ]
            )

        self.energy_bins = np.logspace(2, 3.8, self.spectrum_bins)
        self.energy_out = np.zeros((len(self.energy_bins - 1), self.time_steps))
        self.packets_info_array = np.zeros((int(self.num_packets), 8))

    def time_gamma_packet_loop(self):
        gamma_packet_loop(
            self.packets,
            self.grey_opacity,
            self.photoabsorption_opacity,
            self.pair_creation_opacity,
            self.electron_number_density_time,
            self.mass_density_time,
            self.inv_volume_time,
            self.iron_group_fraction_per_shell,
            self.inner_velocities,
            self.outer_velocities,
            self.times,
            self.dt_array,
            self.effective_time_array,
            self.energy_bins,
            self.energy_df_rows,
            self.energy_plot_df_rows,
            self.energy_out,
            self.packets_info_array
        )
