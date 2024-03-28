"""
Basic TARDIS Benchmark.
"""
import astropy.units as u
from asv_runner.benchmarks.mark import skip_benchmark

from benchmarks.benchmark_base import BenchmarkBase
from tardis.energy_input.energy_source import get_all_isotopes, setup_input_energy
from tardis.energy_input.gamma_ray_transport import create_isotope_dicts, create_inventories_dict, \
    calculate_total_decays, decay_chain_energies, calculate_average_energies, calculate_energy_per_mass


# @skip_benchmark
class BenchmarkEnergyInputGammaRayTransport(BenchmarkBase):
    """
    Class to benchmark the gamma ray transport function.
    """

    def __init__(self):
        pass

    def time_calculate_cell_masses(self):
        """
        Function to test calculation of shell masses.
        """
        self.gamma_ray_simulation_state.composition.calculate_cell_masses(
            self.gamma_ray_simulation_state.geometry.volume
        )

    def time_calculate_total_decays_activity(self):
        # setup of decay test
        time_delta = 1.0 * u.s

        # calculating necessary values
        composition = self.gamma_ray_simulation_state.composition
        cell_masses = composition.calculate_cell_masses(
            self.gamma_ray_simulation_state.geometry.volume
        )
        isotopic_mass_fractions = (
            self.gamma_ray_simulation_state.composition.isotopic_mass_fraction
        )
        iso_dict = create_isotope_dicts(isotopic_mass_fractions, cell_masses)
        inv_dict = create_inventories_dict(iso_dict)

        calculate_total_decays(inv_dict, time_delta)

    def time_calculate_total_decays_activity_chain(self):
        time_delta = 1.0 * u.d.to(u.s)

        composition = self.gamma_ray_simulation_state.composition
        cell_masses = composition.calculate_cell_masses(
            self.gamma_ray_simulation_state.geometry.volume
        )
        isotopic_mass_fractions = (
            self.gamma_ray_simulation_state.composition.isotopic_mass_fraction
        )
        iso_dict = create_isotope_dicts(isotopic_mass_fractions, cell_masses)
        inv_dict = create_inventories_dict(iso_dict)

        calculate_total_decays(inv_dict, time_delta)

    def time_isotope_dicts(self):
        isotopic_mass_fractions = (
            self.gamma_ray_simulation_state.composition.isotopic_mass_fraction
        )
        composition = self.gamma_ray_simulation_state.composition
        cell_masses = composition.calculate_cell_masses(
            self.gamma_ray_simulation_state.geometry.volume
        )
        create_isotope_dicts(isotopic_mass_fractions, cell_masses)

    def time_average_energies(self):
        isotopic_mass_fraction = (
            self.gamma_ray_simulation_state.composition.isotopic_mass_fraction
        )
        gamma_ray_lines = self.atomic_dataset.decay_radiation_data

        all_isotope_names = get_all_isotopes(isotopic_mass_fraction)

        for isotope_name in all_isotope_names:
            setup_input_energy(
                gamma_ray_lines[
                    gamma_ray_lines.index == isotope_name.replace("-", "")
                    ],
                "g",
            )

    def time_decay_energy_chain(self):
        isotopic_mass_fractions = (
            self.gamma_ray_simulation_state.composition.isotopic_mass_fraction
        )

        composition = self.gamma_ray_simulation_state.composition
        cell_masses = composition.calculate_cell_masses(
            self.gamma_ray_simulation_state.geometry.volume
        )
        iso_dict = create_isotope_dicts(isotopic_mass_fractions, cell_masses)
        inventories_dict = create_inventories_dict(iso_dict)
        gamma_ray_lines = self.atomic_dataset.decay_radiation_data

        total_decays = calculate_total_decays(inventories_dict, 1.0 * u.s)

        (
            average_energies,
            _,
            _,
        ) = calculate_average_energies(isotopic_mass_fractions, gamma_ray_lines)

        decay_chain_energies(
            average_energies,
            total_decays,
        )

    def time_energy_per_mass(self):
        isotopic_mass_fractions = (
            self.gamma_ray_simulation_state.composition.isotopic_mass_fraction
        )
        composition = self.gamma_ray_simulation_state.composition
        cell_masses = composition.calculate_cell_masses(
            self.gamma_ray_simulation_state.geometry.volume
        )
        iso_dict = create_isotope_dicts(isotopic_mass_fractions, cell_masses)
        inventories_dict = create_inventories_dict(iso_dict)
        total_decays = calculate_total_decays(inventories_dict, 1.0 * u.s)

        gamma_ray_lines = self.atomic_dataset.decay_radiation_data
        average_energies = calculate_average_energies(
            isotopic_mass_fractions, gamma_ray_lines
        )
        decay_energy = decay_chain_energies(
            average_energies[0],
            total_decays,
        )
        calculate_energy_per_mass(
            decay_energy, isotopic_mass_fractions, cell_masses
        )
