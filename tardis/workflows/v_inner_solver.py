import logging

import numpy as np
from astropy import units as u

from tardis.io.atom_data.base import AtomData
from tardis.simulation.convergence import ConvergenceSolver
from tardis.spectrum.formal_integral import FormalIntegrator
from tardis.workflows.simple_simulation import SimpleSimulation
from tardis.workflows.util import get_tau_integ
from tardis.plasma.standard_plasmas import assemble_plasma
from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.opacities.opacity_solver import OpacitySolver
from tardis.opacities.macro_atom.macroatom_solver import MacroAtomSolver
from tardis.opacities.macro_atom.macroatom_state import MacroAtomState
from scipy.interpolate import interp1d
import copy
from tardis.spectrum.luminosity import (
    calculate_filtered_luminosity,
)
import pandas as pd

# logging support
logger = logging.getLogger(__name__)

# TODO:
# Option to estimate initial v_inner from electron opacity
# Add option for desired optical depth
# Think about csvy vs other formats for specifying grid
# Handle non-explicit formats when going out of the simulation


class InnerVelocitySimulationSolver(SimpleSimulation):

    TAU_TARGET = np.log(2.0 / 3)

    def __init__(
        self,
        configuration,
        mean_optical_depth="rossland",
        tau=None
    ):
        """
        Args:
            convergence_strategy (_type_): _description_
            atom_data_path (_type_): _description_
            mean_optical_depth (str): 'rossland' or 'planck'
                Method of estimating the mean optical depth
        """
        super().__init__(configuration)
        self.mean_optical_depth = mean_optical_depth

        self.convergence_solvers["v_inner_boundary"] = ConvergenceSolver(
            self.convergence_strategy.v_inner_boundary
        )

        self.property_mask_ = self.property_mask

        self.opacity_solver = OpacitySolver(line_interaction_type=configuration.plasma.line_interaction_type, disable_line_scattering=False)
        if configuration.plasma.line_interaction_type == 'scatter':
            self.macro_atom_solver = None
        else:
            self.macro_atom_solver = MacroAtomSolver()
        if tau is not None:
            self.TAU_TARGET = np.log(tau)

    def estimate_v_inner(self):
        """Compute the Rossland Mean Optical Depth,
        Estimate location where v_inner makes t=2/3 (or target)
        Extrapolate with exponential fits

        Need some way to return and inspect the optical depths for later logging"""
        pass

        tau_integ = np.log(
            get_tau_integ(
                self.plasma_solver,
                self.simulation_state,
            )[self.mean_optical_depth]
        )

        interpolator = interp1d(
            tau_integ,
            self.simulation_state.geometry.v_inner,  # Only use the active values as we only need a numerical estimate, not an index
            fill_value="extrapolate",
        )
        # TODO: Make sure eastimed_v_inner is within the bounds of the simulation!
        estimated_v_inner = interpolator(self.TAU_TARGET)

        return estimated_v_inner * u.cm / u.s

    @property
    def property_mask(self):
        mask = np.zeros(
            (len(self.simulation_state.geometry.r_inner)), dtype=bool
        )
        mask[
            self.simulation_state.geometry.v_inner_boundary_index : self.simulation_state.geometry.v_outer_boundary_index
        ] = True
        return mask

    @property
    def property_where(self):

        return np.where(self.property_mask)

    def get_convergence_estimates(self, transport_state):

        # print(self.transport_solver.transport_state)
        # (
        #    estimated_t_radiative,
        #    estimated_dilution_factor,
        # ) = self.transport_solver.radfield_prop_solver.calculate_radiationfield_properties()
        
        estimated_radfield_properties = (
            self.transport_solver.radfield_prop_solver.solve(
                transport_state.radfield_mc_estimators,
                transport_state.time_explosion,
                transport_state.time_of_simulation,
                transport_state.geometry_state.volume,
                transport_state.opacity_state.line_list_nu,
            )
        )
        print("Volume_INNER:", transport_state.geometry_state.volume[0])

        estimated_t_radiative = estimated_radfield_properties.dilute_blackbody_radiationfield_state.temperature
        estimated_dilution_factor = estimated_radfield_properties.dilute_blackbody_radiationfield_state.dilution_factor
        
        self.initialize_spectrum_solver(
            transport_state,
            None,
        )

        emitted_luminosity = calculate_filtered_luminosity(
            transport_state.emitted_packet_nu,
            transport_state.emitted_packet_luminosity,
            self.luminosity_nu_start,
            self.luminosity_nu_end,
        )

        print("Emitted Luminosity:", emitted_luminosity)
        luminosity_ratios = (
            (emitted_luminosity / self.luminosity_requested).to(1).value
        )

        estimated_t_inner = (
            self.simulation_state.t_inner
            * luminosity_ratios
            ** self.convergence_strategy.t_inner_update_exponent
        )

        self.tracker['emitted_luminosity'].append(emitted_luminosity)
        self.tracker['estimated_t_inner'].append(estimated_t_inner)

        estimated_v_inner = self.estimate_v_inner()
        if estimated_v_inner < self.simulation_state.geometry.v_inner[0]:
            estimated_v_inner = self.simulation_state.geometry.v_inner[0]
            print("WARNING: v_inner_boundary outside of simulation, setting to first shell")
        elif estimated_v_inner > self.simulation_state.geometry.v_inner[-1]:
            estimated_v_inner = self.simulation_state.geometry.v_inner[-1]
            print("WARNING: v_inner_boundary outside of simulation, setting to last shell")


        #estimated_v_inner = self.simulation_state.v_inner_boundary
        print(estimated_v_inner)

        return {
            "t_radiative": estimated_t_radiative,
            "dilution_factor": estimated_dilution_factor,
            "t_inner": estimated_t_inner,
            "v_inner_boundary": estimated_v_inner,
        }
    
    def reproject(self, a1, a2, m1, m2):

        a1_expanded = np.empty(
        len(self.simulation_state.geometry.r_inner),
            dtype=a1.dtype,
        )
        a2_expanded = np.empty_like(a1_expanded)

        a1_expanded[m1] = a1
        a2_expanded[m2] = a2
        
        joint_mask = m1 & m2

        return a1_expanded[joint_mask], a2_expanded[joint_mask]


    def check_convergence(
        self,
        estimated_values,
    ):
        convergence_statuses = []

        for key, solver in self.convergence_solvers.items():

            current_value = getattr(self.simulation_state, key)
            estimated_value = estimated_values[key]

            print("Check Convergence")
            print(key, estimated_value)
            print(current_value)
            joint_mask = np.ones(5)
            if hasattr(current_value, "__len__") and (
                key not in ["t_inner", "v_inner_boundary"]
            ):
                if current_value.shape == estimated_value.shape:
                    pass
                else:
                    print(key, "has", "__len__")
                    print("shape current:", current_value.shape)
                    print("shape next:", estimated_value.shape)
                    new_value = estimated_value
                    current_value_expanded = np.empty(
                        len(self.simulation_state.geometry.r_inner),
                        dtype=current_value.dtype,
                    )
                    current_value_expanded[self.property_mask] = current_value
                    new_value_expanded = np.empty_like(current_value_expanded)
                    new_value_expanded[self.new_property_mask] = new_value
                    joint_mask = self.property_mask & self.new_property_mask
                    print(joint_mask)
                    print("diff:", current_value_expanded - new_value_expanded)

                    if hasattr(current_value, "unit"):
                        current_value_expanded = (
                            current_value_expanded * current_value.unit
                        )
                        new_value_expanded = (
                            new_value_expanded * current_value.unit
                        )
                    estimated_value = new_value_expanded[joint_mask]
                    current_value = current_value_expanded[joint_mask]

            # no_of_shells = (
            #    self.simulation_state.no_of_shells if key not in ["t_inner", "v_inner_boundary"] else 1
            # )
            no_of_shells = (
                joint_mask.sum()
                if key not in ["t_inner", "v_inner_boundary"]
                else 1
            )

            convergence_statuses.append(
                solver.get_convergence_status(
                    current_value, estimated_value, no_of_shells
                )
            )
            print("Status:", convergence_statuses[-1])

        if np.all(convergence_statuses):
            hold_iterations = self.convergence_strategy.hold_iterations
            self.consecutive_converges_count += 1
            logger.info(
                f"Iteration converged {self.consecutive_converges_count:d}/{(hold_iterations + 1):d} consecutive "
                f"times."
            )
            print("Converged this iteration!")
            return self.consecutive_converges_count == hold_iterations + 1

        self.consecutive_converges_count = 0
        return False

    def clip(self, property):
        """Clips a shell-dependent array to the current index"""

        return property[
            self.simulation_state.geometry.v_inner_boundary_index : self.simulation_state.geometry.v_outer_boundary_index
        ]

    def solve_simulation_state(
        self,
        estimated_values,
    ):

        next_values = {}
        print(estimated_values)
        self.new_property_mask = self.property_mask
        self.old_property_mask = self.property_mask_

        for key, solver in self.convergence_solvers.items():
            if (
                key in ["t_inner", "v_inner_boundary"]
                and (self.completed_iterations + 1)
                % self.convergence_strategy.lock_t_inner_cycles
                != 0
            ):
                next_values[key] = getattr(self.simulation_state, key)
            else:

                print("key", key)
                print(getattr(self.simulation_state, key))
                print(estimated_values[key])
                current_value = getattr(self.simulation_state, key)
                new_value = estimated_values[key]
                if hasattr(current_value, "__len__") and key not in [
                    "t_inner",
                    "v_inner_boundary",
                ]:
                    if current_value.shape == new_value.shape:
                        pass
                    else:
                        breakpoint()
                        print(key, "has", "__len__")
                        print("shape current:", current_value.shape)
                        print("shape next:", new_value.shape)
                        current_value_expanded = np.empty(
                            len(self.simulation_state.geometry.r_inner),
                            dtype=current_value.dtype,
                        )
                        current_value_expanded[
                            self.property_mask
                        ] = current_value
                        new_value_expanded = np.empty_like(
                            current_value_expanded
                        )
                        new_value_expanded[self.new_property_mask] = new_value
                        joint_mask = self.property_mask & self.new_property_mask
                        print(joint_mask)
                        if hasattr(current_value, "unit"):
                            current_value_expanded = (
                                current_value_expanded * current_value.unit
                            )
                            new_value_expanded = (
                                new_value_expanded * current_value.unit
                            )
                        new_value = new_value_expanded[joint_mask]
                        current_value = current_value_expanded[joint_mask]
                next_values[key] = solver.converge(
                    current_value, new_value
                )  # TODO: This needs to be changed to account for changing array sizes

        self.simulation_state.t_radiative = next_values["t_radiative"]
        self.simulation_state.dilution_factor = next_values["dilution_factor"]
        self.simulation_state.blackbody_packet_source.temperature = next_values[
            "t_inner"
        ]
        self.simulation_state.t_inner = next_values["t_inner"]

        print("next v_inner", next_values["v_inner_boundary"])
        self.simulation_state.geometry.v_inner_boundary = next_values[
            "v_inner_boundary"
        ]
        self.simulation_state.blackbody_packet_source.radius = self.simulation_state.r_inner[0]
        self.property_mask_ = self.new_property_mask

        print("New Volume:", self.simulation_state.geometry.volume_active)


    def solve_montecarlo(self, no_of_real_packets, no_of_virtual_packets=0):

        opacity_state = self.opacity_solver.solve(self.plasma_solver)
        #macro_atom_state = None
        if self.macro_atom_solver is None:
            macro_atom_state = None
        else:
            macro_atom_state = self.macro_atom_solver.solve(
                    self.plasma_solver,
                    self.plasma_solver.atomic_data,
                    opacity_state.tau_sobolev,
                    self.plasma_solver.stimulated_emission_factor,
                )
        transport_state = self.transport_solver.initialize_transport_state(
            self.simulation_state,
            opacity_state,
            macro_atom_state,
            self.plasma_solver,
            no_of_real_packets,
            no_of_virtual_packets=no_of_virtual_packets,
            iteration=self.completed_iterations,
        )

        virtual_packet_energies = self.transport_solver.run(
            transport_state,
            iteration=self.completed_iterations,
            total_iterations=self.total_iterations,
            show_progress_bars=False,
        )

        output_energy = transport_state.packet_collection.output_energies
        if np.sum(output_energy < 0) == len(output_energy):
            logger.critical("No r-packet escaped through the outer boundary.")

        return transport_state, virtual_packet_energies

    def solve_plasma(
        self,
        transport_state,
    ):
        # TODO: Find properties that need updating with shells

        # self.simulation_state.radiation_field_state.t_radiative[self.simulation_state.radiation_field_state.t_radiative.value==0] = 10000.0 * u.K
        # self.simulation_state.radiation_field_state.dilution_factor[self.simulation_state.radiation_field_state.dilution_factor==0] = 1.0

        radiation_field = DilutePlanckianRadiationField(
            temperature=self.simulation_state.radiation_field_state.temperature,
            dilution_factor=self.simulation_state.radiation_field_state.dilution_factor,
        )
        update_properties = dict(
            dilute_planckian_radiation_field=radiation_field,
        )
        # A check to see if the plasma is set with JBluesDetailed, in which
        # case it needs some extra kwargs.
        j_blues = radiation_field.calculate_mean_intensity(
                self.plasma_solver.atomic_data.lines.nu.values
            )
        update_properties["j_blues"] = pd.DataFrame(
            j_blues, index=self.plasma_solver.atomic_data.lines.index
        )
        if "j_blue_estimator" in self.plasma_solver.outputs_dict:
            update_properties.update(
                t_inner=self.simulation_state.blackbody_packet_source.temperature,
                j_blue_estimator=transport_state.radfield_mc_estimators.j_blue_estimator,
            )

        self.plasma_solver.update(**update_properties)

    def solve(self):
        converged = False
        while self.completed_iterations < self.total_iterations - 1:

            

            transport_state, virtual_packet_energies = self.solve_montecarlo(
                self.real_packet_count
            )

            estimated_values = self.get_convergence_estimates(transport_state)

            self.solve_simulation_state(estimated_values)

            self.solve_plasma(transport_state)

            converged = self.check_convergence(estimated_values)

            #self.simulation_state.packet_source.radius = self.simulation_state.geometry.r_inner_active[0]

            self.completed_iterations += 1
            if converged:
                print("SIMULATION CONVERGED!")
            if converged and self.convergence_strategy.stop_if_converged:
                break

        transport_state, virtual_packet_energies = self.solve_montecarlo(
            self.final_iteration_packet_count, self.virtual_packet_count
        )
        self.initialize_spectrum_solver(
            transport_state,
            virtual_packet_energies,
        )
