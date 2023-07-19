import radioactivedecay as rd
from radioactivedecay.utils import Z_DICT
from tardis.io.config_reader import Configuration
from tardis.model.base import ModelState, Composition
from tardis.model.geometry.radial1d import Radial1DGeometry
from tardis.radiation_field.source_functions.base import PhotosphericBlackBody1D
from tardis.montecarlo.packet_source import (
    BlackBodySimpleSource,
    BlackBodySimpleSourceRelativistic,
)
from tardis.montecarlo.base import MontecarloTransport
from tardis.radiation_field.base import MonteCarloRadiationFieldState
from tardis.io.atom_data.base import AtomData

from tardis.model.density import (
    HomologousDensityState,
    calculate_power_law_density,
    calculate_exponential_density,
)
from tardis.util.base import quantity_linspace
import pandas as pd
from tardis.io.model_reader import (
    read_density_file,
    read_abundances_file,
    read_uniform_abundances,
    parse_csv_abundances,
)
import os
import numpy as np
from astropy import units as u
from astropy import constants as const
import logging
from tardis.io.decay import IsotopeAbundances
from tardis.simulation.base import Simulation
from tardis.plasma.standard_plasmas import assemble_plasma

logger = logging.getLogger(__name__)

GRAPH_EXPLORE_CALLS = 0 # For pure curiosity
def requires_single(attr: str):
    def decorator(func):
        def _wrapper(self, *args, **kwargs):
            global GRAPH_EXPLORE_CALLS
            if not hasattr(self, attr):
                GRAPH_EXPLORE_CALLS += 1
                if not hasattr(self, f"to_{attr}"):
                    logger.warn(
                        f"No attribute {attr} provided, setting {attr}=None"
                    )
                    self.__setitem__(attr, None)
                else:
                    exec(f"self.to_{attr}()")
            return func(self, *args, **kwargs)

        # Make sure our returned function looks just like our original function
        _wrapper.__doc__ = func.__doc__
        if hasattr(func, "__qualname__"):
            _wrapper.__qualname__ = func.__qualname__
        return _wrapper

    return decorator


def requires(*attrs):
    def decorator(func):
        new_func = func
        for attr in attrs:
            new_func = requires_single(attr)(new_func)
        return new_func

    return decorator


# This can either be subclasses, reclassed, or extra ways to parse like with csvy can be added externally to the from_config_dict
# Only currently works for model.structure.type=='specific'
class TARDISFactory(Configuration):


    @classmethod
    def from_config_dict(
        cls, config_dict, validate=True, config_dirname="", **kwargs
    ):
        new_config = super().from_config_dict(
            config_dict, validate=validate, config_dirname=config_dirname
        )
        new_config.__dict__.update(kwargs)
        return new_config

    def to_velocity(self):
        
        self.velocity = quantity_linspace(
            self.model.structure.velocity.start,
            self.model.structure.velocity.stop,
            self.model.structure.velocity.num + 1,
        ).cgs

        adjusted_velocity = self.velocity.insert(0, 0)
        self.v_middle = (
            adjusted_velocity[1:] * 0.5 + adjusted_velocity[:-1] * 0.5
        )
        return self.velocity

    @requires("velocity")
    def to_branch85_w7(self):
        density_0 = calculate_power_law_density(
            self.v_middle,
            self.model.structure.density.w7_v_0,
            self.model.structure.density.w7_rho_0,
            -7,
        )
        time_0 = self.model.structure.density.w7_time_0
        self.branch85_w7 = HomologousDensityState(density_0, time_0)
        return self.branch85_w7

    @requires("velocity")
    def to_exponential(self):
        density_0 = calculate_exponential_density(
            self.v_middle,
            self.model.structure.density.v_0,
            self.model.structure.density.rho_0,
        )
        time_0 = self.model.structure.density.time_explosion
        self.exponential = HomologousDensityState(density_0, time_0)
        return self.exponential

    @requires("velocity")
    def to_power_law(self):
        density_0 = calculate_power_law_density(
            self.v_middle,
            self.model.structure.density.v_0,
            self.model.structure.density.rho_0,
            self.model.structure.density.exponent,
        )
        time_0 = self.model.structure.density.get(
            "time_0", self.model.structure.density.time_explosion
        )
        self.power_law = HomologousDensityState(density_0, time_0)
        return self.power_law

    @requires("velocity", "no_of_shells")
    def to_uniform(self):
        density_0 = self.model.structure.density.value.to("g cm^-3") * np.ones(
            self.no_of_shells
        )
        time_0 = self.model.structure.density.get(
            "time_0", self.model.structure.density.time_explosion
        )
        self.uniform = HomologousDensityState(density_0, time_0)
        return self.uniform

    @requires("geometry")
    def to_no_of_shells(self):
        self.no_of_shells = len(self.geometry.v_inner)
        return self.no_of_shells

    @requires(
        "atom_data"
    )  # This method and associated property is poorly worded
    def to_atomic_data(self):
        self.atomic_data = AtomData.from_hdf(self.atom_data)
        return self.atomic_data

    @requires("model_state", "atomic_data", "radiation_field")
    def to_plasma_state(self):
        # This is going to be a little cheat here as a test, needs more refactoring
        self.model_state._electron_densities = None
        self.plasma_state = assemble_plasma(
            self,
            self.model_state,
            self.radiation_field,
            atom_data=self.atomic_data,
        )
        return self.plasma_state

    @requires(
        "model_state", "radiation_field", "transport_state", "plasma_state"
    )
    def to_simulation(
        self, show_convergence_plots=True, show_progress_bars=True, **kwargs
    ):
        luminosity_nu_start = self.supernova.luminosity_wavelength_end.to(
            u.Hz, u.spectral()
        )

        if u.isclose(
            self.supernova.luminosity_wavelength_start, 0 * u.angstrom
        ):
            luminosity_nu_end = np.inf * u.Hz
        else:
            luminosity_nu_end = (
                const.c / self.supernova.luminosity_wavelength_start
            ).to(u.Hz)

        last_no_of_packets = self.montecarlo.last_no_of_packets
        if last_no_of_packets is None or last_no_of_packets < 0:
            last_no_of_packets = self.montecarlo.no_of_packets
        last_no_of_packets = int(last_no_of_packets)

        convergence_plots_config_options = [
            "plasma_plot_config",
            "t_inner_luminosities_config",
            "plasma_cmap",
            "t_inner_luminosities_colors",
            "export_convergence_plots",
        ]
        convergence_plots_kwargs = {}
        for item in set(convergence_plots_config_options).intersection(
            kwargs.keys()
        ):
            convergence_plots_kwargs[item] = kwargs[item]

        Simulation(
            iterations=self.montecarlo.iterations,
            model_state=self.model_state,
            plasma=self.plasma_state,
            radiation_field=self.radiation_field,
            transport=self.transport_state,
            show_convergence_plots=show_convergence_plots,
            no_of_packets=int(self.montecarlo.no_of_packets),
            no_of_virtual_packets=int(self.montecarlo.no_of_virtual_packets),
            luminosity_nu_start=luminosity_nu_start,
            luminosity_nu_end=luminosity_nu_end,
            last_no_of_packets=last_no_of_packets,
            luminosity_requested=self.supernova.luminosity_requested.cgs,
            convergence_strategy=self.montecarlo.convergence_strategy,
            convergence_plots_kwargs=convergence_plots_kwargs,
            show_progress_bars=show_progress_bars,
        )

    @requires("geometry", "composition", "time_explosion")
    def to_model_state(self):
        self.model_state = ModelState(
            self.composition, self.geometry, self.time_explosion
        )
        return self.model_state

    def to_virtual_packet_logging(self):
        self.virtual_packet_logging = (
            self.spectrum.virtual.virtual_packet_logging
        )
        return self.virtual_packet_logging

    @requires("source_function", "virtual_packet_logging")
    def to_transport_state(self):
        self.transport_state = MontecarloTransport.from_config(
            config,
            packet_source=self.source_function.packet_source,
            virtual_packet_logging=self.virtual_packet_logging,
        )
        return self.transport_state

    @requires("density", "atomic_data", "no_of_shells")
    def to_composition(self):
        abundances_section = self.model.abundances
        isotope_abundance = pd.DataFrame()
        if abundances_section.type == "uniform":  # The only supported scheme
            (
                elemental_mass_fraction,
                isotope_abundance,
            ) = read_uniform_abundances(abundances_section, self.no_of_shells)

        else:
            raise Exception("Only Uniform Abundances Supported")

        elemental_mass_fraction = elemental_mass_fraction.replace(np.nan, 0.0)
        elemental_mass_fraction = elemental_mass_fraction[
            elemental_mass_fraction.sum(axis=1) > 0
        ]

        norm_factor = elemental_mass_fraction.sum(
            axis=0
        ) + isotope_abundance.sum(axis=0)

        if np.any(np.abs(norm_factor - 1) > 1e-12):
            logger.warning(
                "Abundances have not been normalized to 1." " - normalizing"
            )
            elemental_mass_fraction /= norm_factor
            isotope_abundance /= norm_factor

        isotope_abundance = IsotopeAbundances(isotope_abundance)

        atomic_mass = None
        elemental_mass = self.atomic_data.atom_data.mass

        mass = {}
        stable_atomic_numbers = elemental_mass_fraction.index.to_list()

        for z in stable_atomic_numbers:
            mass[z] = [
                elemental_mass[z]
                for i in range(elemental_mass_fraction.columns.size)
            ]
        stable_isotope_mass = pd.DataFrame(mass).T

        isotope_mass = {}
        for atomic_number, i in isotope_abundance.decay(
            self.time_explosion
        ).groupby(level=0):
            i = i.loc[atomic_number]
            for column in i:
                mass = {}
                shell_mass_fractions = i[column]
                isotopic_masses = [
                    rd.Nuclide(Z_DICT[atomic_number] + str(i)).atomic_mass
                    for i in shell_mass_fractions.index.to_numpy()
                ]
                mass[atomic_number] = (
                    shell_mass_fractions * isotopic_masses
                ).sum()
                mass[atomic_number] /= shell_mass_fractions.sum()
                mass[atomic_number] = mass[atomic_number] * u.u.to(u.g)
                if isotope_mass.get(column) is None:
                    isotope_mass[column] = {}
                isotope_mass[column][atomic_number] = mass[atomic_number]
        isotope_mass = pd.DataFrame(isotope_mass)

        atomic_mass = pd.concat([stable_isotope_mass, isotope_mass])

        self.composition = Composition(
            self.density,
            elemental_mass_fraction,
            atomic_mass,
            atomic_mass_unit=u.g,
        )
        return self.composition

    @requires("source_function", "transport_state", "no_of_shells")
    def to_radiation_field(self):
        # Using's Wolfgang's clever estimation scheme
        lambda_wien_inner = const.b_wien / self.source_function.t_inner
        t_radiative = const.b_wien / (
            lambda_wien_inner
            * (1 + (self.v_middle - self.geometry.v_inner[0]) / const.c)
        )
        r_middle = (self.geometry.r_inner + self.geometry.r_outer) / 2
        dilution_factor = 0.5 * (
            1
            - np.sqrt(
                1 - (self.geometry.r_inner[0] ** 2 / r_middle**2).to(1).value
            )
        )
        opacities = None
        self.radiation_field = MonteCarloRadiationFieldState(
            t_radiative, dilution_factor, opacities, self.source_function
        )
        return self.radiation_field

    @requires("geometry")
    def to_source_function(self):
        v_inner_boundary = self.geometry.v_inner[0]
        t_inner = self.plasma.initial_t_inner
        no_of_packets = self.montecarlo.no_of_packets

        if not hasattr(self, "packet_source") is None:
            if not self.montecarlo.enable_full_relativity:
                self.packet_source = BlackBodySimpleSource(self.montecarlo.seed)
            else:
                self.packet_source = BlackBodySimpleSourceRelativistic(
                    self.montecarlo.seed
                )

        self.source_function = PhotosphericBlackBody1D(
            v_inner_boundary, t_inner, no_of_packets, self.packet_source
        )
        return self.source_function

    def to_density(self):
        """
        Create a new HomologousDensityState instance from a Configuration object.

        Parameters
        ----------
        config : tardis.io.config_reader.Configuration

        Returns
        -------
        HomologousDensityState

        """

        # This also serves as an example of how to handle multiple options for how to parse something
        # Could also have the .to_density method map the parsing method to the parser as an alternative
        @requires(self.model.structure.density.type)
        def get_density(self):
            self.density = self.__getattr__(self.model.structure.density.type)
            return self.density

        return get_density(self)

    def to_time_explosion(self):
        self.time_explosion = self.supernova.time_explosion.cgs
        return self.time_explosion

    @requires("velocity", "time_explosion")
    def to_geometry(self):
        v_outer = self.velocity[1:]
        v_inner = self.velocity[:-1]
        r_inner = v_inner * self.time_explosion
        r_outer = v_outer * self.time_explosion

        self.geometry = Radial1DGeometry(r_inner, r_outer, v_inner, v_outer)

        return self.geometry

    def to_tardis(self):
        pass


if __name__ == "__main__":
    # For testing
    config = TARDISFactory.from_yaml("../../docs/tardis_example.yml")
    density = config.to_density()
    geometry = config.to_geometry()
    sim = config.to_simulation() # NEED TO MAKE SURE TO EVOLVE THE MODEL STATE!
    print(GRAPH_EXPLORE_CALLS) 
    breakpoint()
