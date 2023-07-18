from tardis.io.config_reader import Configuration, q
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
from tardis.plasma.base import assemble_plasma
logger = logging.getLogger(__name__)

class TARDISFactory(Configuration):


    def to_simulation(self, 
        packet_source=None,
        virtual_packet_logging=False,
        show_convergence_plots=True,
        show_progress_bars=True,
        **kwargs):

        if "model" in kwargs:
            model_state = kwargs["model"]
        else:
            model_state = self.to_model_state()
            radiation_field = self.to_radiation_field()
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

        pass

    def to_transport(self):

        pass

    def to_composition(self):

       pass

    def to_radiation_field(self):

        pass

    def to_packet_source(self, packet_source=None):


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
        

    def to_density_csvy(self):

        """
        Create a new HomologousDensity instance from a base
        Configuration object and a csvy model Configuration object.

        Parameters
        ----------
        config : tardis.io.config_reader.Configuration
        csvy_model_config : tardis.io.config_reader.Configuration

        Returns
        -------
        HomologousDensity

        """
        if hasattr(csvy_model_config, "velocity"):
            velocity = quantity_linspace(
                csvy_model_config.velocity.start,
                csvy_model_config.velocity.stop,
                csvy_model_config.velocity.num + 1,
            ).cgs
        else:
            velocity_field_index = [
                field.name for field in csvy_model_config.datatype.fields
            ].index("velocity")
            velocity_unit = u.Unit(
                csvy_model_config.datatype.fields[velocity_field_index].unit
            )
            velocity = csvy_model_config.velocity.values * velocity_unit

        adjusted_velocity = velocity.insert(0, 0)
        v_middle = adjusted_velocity[1:] * 0.5 + adjusted_velocity[:-1] * 0.5
        no_of_shells = len(adjusted_velocity) - 1
        time_explosion = config.supernova.time_explosion.cgs

        if hasattr(csvy_model_config, "density"):
            d_conf = csvy_model_config.density
            density_type = d_conf.type
            if density_type == "branch85_w7":
                density_0 = calculate_power_law_density(
                    v_middle, d_conf.w7_v_0, d_conf.w7_rho_0, -7
                )
                time_0 = d_conf.w7_time_0
            elif density_type == "uniform":
                density_0 = d_conf.value.to("g cm^-3") * np.ones(no_of_shells)
                time_0 = d_conf.get("time_0", time_explosion)
            elif density_type == "power_law":
                density_0 = calculate_power_law_density(
                    v_middle, d_conf.v_0, d_conf.rho_0, d_conf.exponent
                )
                time_0 = d_conf.get("time_0", time_explosion)
            elif density_type == "exponential":
                density_0 = calculate_exponential_density(
                    v_middle, d_conf.v_0, d_conf.rho_0
                )
                time_0 = d_conf.get("time_0", time_explosion)
            else:
                raise ValueError(
                    f"Unrecognized density type " f"'{d_conf.type}'"
                )
        return cls(density_0, time_0)

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
        
        d_conf = self.model.structure.density
        velocity = quantity_linspace(
            self.model.structure.velocity.start,
            self.model.structure.velocity.stop,
            self.model.structure.velocity.num + 1,
        ).cgs

        adjusted_velocity = velocity.insert(0, 0)
        v_middle = adjusted_velocity[1:] * 0.5 + adjusted_velocity[:-1] * 0.5
        no_of_shells = len(adjusted_velocity) - 1
        time_explosion = self.supernova.time_explosion.cgs

        if d_conf.type == "branch85_w7":
            density_0 = calculate_power_law_density(
                v_middle, d_conf.w7_v_0, d_conf.w7_rho_0, -7
            )
            time_0 = d_conf.w7_time_0
        elif d_conf.type == "uniform":
            density_0 = d_conf.value.to("g cm^-3") * np.ones(no_of_shells)
            time_0 = d_conf.get("time_0", time_explosion)
        elif d_conf.type == "power_law":
            density_0 = calculate_power_law_density(
                v_middle, d_conf.v_0, d_conf.rho_0, d_conf.exponent
            )
            time_0 = d_conf.get("time_0", time_explosion)
        elif d_conf.type == "exponential":
            density_0 = calculate_exponential_density(
                v_middle, d_conf.v_0, d_conf.rho_0
            )
            time_0 = d_conf.get("time_0", time_explosion)
        else:
            raise ValueError(f"Unrecognized density type " f"'{d_conf.type}'")
        return HomologousDensityState(density_0, time_0)

    def to_geometry(self):


        pass

    def to_tardis(self):

        pass