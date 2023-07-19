import radioactivedecay as rd
from radioactivedecay.utils import Z_DICT
from tardis.io.config_reader import Configuration
from tardis.model.base import ModelState, Composition
from tardis.model.geometry.radial1d import Radial1DGeometry
from tardis.radiation_field.source_functions.base import PhotosphericBlackBody1D
from tardis.montecarlo.packet_source import BlackBodySimpleSource, BlackBodySimpleSourceRelativistic
from tardis.model.density import HomologousDensityState, calculate_power_law_density, calculate_exponential_density
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

def requires(attr):

    def decorator(func):

        def _wrapper(self, *args, **kwargs):
            if hasattr(self, attr):
                pass
            else:
                exec(f'self.to_{attr}()')
            return func(self, *args, **kwargs)
        return _wrapper
    
    return decorator


class TARDISFactory(Configuration):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.density_parsers = {}

    def to_velocity(self):

        self.velocity = quantity_linspace(
            self.model.structure.velocity.start,
            self.model.structure.velocity.stop,
            self.model.structure.velocity.num + 1,
        ).cgs

        adjusted_velocity = self.velocity.insert(0, 0)
        self.v_middle = adjusted_velocity[1:] * 0.5 + adjusted_velocity[:-1] * 0.5
        return self.velocity

    @requires('velocity')
    def to_branch85_w7(self):
        density_0 = calculate_power_law_density(
                self.v_middle, self.model.structure.density.w7_v_0, self.model.structure.density.w7_rho_0, -7
            )
        time_0 = self.model.structure.density.w7_time_0
        self.branch85_w7 = HomologousDensityState(density_0, time_0)
        return self.branch85_w7

    @requires('velocity')
    def to_exponential(self):
        density_0 = calculate_exponential_density(
                self.v_middle, self.model.structure.density.v_0, self.model.structure.density.rho_0
            )
        time_0 = self.model.structure.density.time_explosion
        self.exponential = HomologousDensityState(density_0, time_0)
        return self.exponential

    @requires('velocity')
    def to_power_law(self):
        density_0 = calculate_power_law_density(
                self.v_middle, self.model.structure.density.v_0, self.model.structure.density.rho_0, self.model.structure.density.exponent
            )
        time_0 = self.model.structure.density.get("time_0", self.model.structure.densitytime_explosion)
        self.power_law = HomologousDensityState(density_0, time_0)
        return self.power_law
        
    @requires('no_of_shells')
    @requires('velocity')
    def to_uniform(self):

        density_0 = self.model.structure.density.value.to("g cm^-3") * np.ones(self.no_of_shells)
        time_0 = self.model.structure.density.get("time_0", self.model.structure.densitytime_explosion)
        self.uniform = HomologousDensityState(density_0, time_0)
        return self.uniform

    @requires('geometry')
    def to_no_of_shells(self):
        self.no_of_shells = len(self.geometry.v_inner)
        return self.no_of_shells

    def to_simulation(self, 
        packet_source=None,
        virtual_packet_logging=False,
        show_convergence_plots=True,
        show_progress_bars=True,
        **kwargs):

        if "model" in kwargs:
            model_state = kwargs["model"]
        else:
            radiation_field = self.to_radiation_field()
            model_state = self.to_model_state()
            
        if "plasma" in kwargs:
            plasma = kwargs["plasma"]
        else:
            plasma = assemble_plasma(
                self, model_state, atom_data=kwargs.get("atom_data", None)
            )
        if "transport" in kwargs:
            transport = kwargs["transport"]
        else:
            transport = self.to_transport()

        if radiation_field.packet_source.v_inner_boundary is None:
            radiation_field.packet_source.v_inner_boundary = model_state.geometry.v_inner[0]


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
            model_state=model_state,
            plasma=plasma,
            radiation_field=radiation_field,
            transport=transport,
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

    def to_model_state(self, atom_data=None):
        
        self.csvy_flag = False
        if hasattr(self, 'csvy_model'):
            self.csvy_flag = True

        time_explosion = self.supernova.time_explosion.cgs
        geometry = self.to_geometry()
        composition = self.to_composition()
        model_state = ModelState(composition, geometry, time_explosion)
        return model_state
        

    def to_transport(self):

        pass

    def to_composition(self):

        if self.csvy_flag:
            density = self.to_density_csvy()
        else:
            density = self.to_density()

        elemental_mass_fraction = None
        atomic_mass = None

        if elemental_mass_fraction is not None:
            mass = {}
            stable_atomic_numbers = self.elemental_mass_fraction.index.to_list()
            for z in stable_atomic_numbers:
                mass[z] = [
                    elemental_mass_fraction[z]
                    for i in range(self.elemental_mass_fraction.columns.size)
                ]
            stable_isotope_mass = pd.DataFrame(mass).T

            isotope_mass = {}
            for atomic_number, i in self.raw_isotope_mass_fraction.decay(
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

        return Composition(density, elemental_mass_fraction, atomic_mass, atomic_mass_unit=u.g)

    def to_radiation_field(self):

        pass

    def to_packet_source(self, packet_source=None, v_boundary_inner=None):


        if hasattr(self, "model"):
            if hasattr(self.model, "v_inner_boundary"):
                v_inner_boundary = self.model.v_inner_boundary
            else:
                v_inner_boundary = None
        else:
            v_inner_boundary = self.model.structure.v_inner_boundary
        t_inner = self.plasma.initial_t_inner
        no_of_packets = self.montecarlo.no_of_packets

        if packet_source is None:
            if not self.montecarlo.enable_full_relativity:
                packet_source = BlackBodySimpleSource(
                    self.montecarlo.seed
                )
            else:
                packet_source = BlackBodySimpleSourceRelativistic(
                    self.montecarlo.seed
                )

        packet_source = PhotosphericBlackBody1D(v_inner_boundary, t_inner, no_of_packets, packet_source=None)
        return packet_source
        
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
        @requires(self.model.structure.density.type)
        def get_density(self):
            self.density = self.__getattr__(self.model.structure.density.type)
            return self.density
        return get_density(self)
    
    def to_time_explosion(self):
        self.time_explosion = self.supernova.time_explosion.cgs
        return self.time_explosion

    @requires('time_explosion')
    @requires('velocity')
    def to_geometry(self):
                
        v_outer = self.velocity[1:]
        v_inner = self.velocity[:-1]
        r_inner = v_inner * self.time_explosion
        r_outer = v_outer * self.time_explosion

        self.geometry = Radial1DGeometry(r_inner, r_outer, v_inner, v_outer)

        return self.geometry

    def to_tardis(self):

        pass

if __name__ == '__main__':
    # For testing
    config = TARDISFactory.from_yaml('../../docs/tardis_example.yml')
    density = config.to_density()
    geometry = config.to_geometry()
    breakpoint()